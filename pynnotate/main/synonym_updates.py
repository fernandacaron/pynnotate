"""Versioned synonym database loading and manual update support."""

import hashlib
import json
import os
import ssl
import sys
import tempfile
import unicodedata
from pathlib import Path
from urllib.parse import urljoin, urlparse
from urllib.request import Request, urlopen


SCHEMA_VERSION = 1
SUPPORTED_DATABASES = ("animal_mito", "animal_nuclear", "plant_mito", "plant_chloro")
DEFAULT_MANIFEST_URL = (
    "https://raw.githubusercontent.com/fernandacaron/pynnotate/"
    "main/pynnotate/data/synonyms_manifest.json"
)
MAX_MANIFEST_BYTES = 64 * 1024
MAX_DATABASE_BYTES = 5 * 1024 * 1024


class SynonymDatabaseError(ValueError):
    """Raised when a synonym database or update is invalid."""


def _normalize_alias(value):
    return unicodedata.normalize("NFKC", value).strip().upper()


def _bundled_database_path():
    candidates = [
        Path(__file__).resolve().parents[1] / "data" / "synonyms_database.json",
    ]
    if getattr(sys, "_MEIPASS", None):
        bundle_root = Path(sys._MEIPASS)
        candidates.extend([
            bundle_root / "data" / "synonyms_database.json",
            bundle_root / "pynnotate" / "data" / "synonyms_database.json",
        ])
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise SynonymDatabaseError("Bundled synonym database was not found.")


def get_synonym_data_dir():
    override = os.environ.get("PYNNOTATE_DATA_DIR", "").strip()
    if override:
        return Path(override).expanduser()
    if sys.platform == "win32":
        base = Path(os.environ.get("APPDATA") or Path.home() / "AppData" / "Roaming")
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support"
    else:
        base = Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share")
    return base / "Pynnotate"


def get_cached_database_path():
    return get_synonym_data_dir() / "synonyms_database.json"


def _validate_remote_url(url):
    parsed = urlparse(url)
    if parsed.scheme == "https" and parsed.netloc:
        return
    if parsed.scheme == "http" and parsed.hostname in {"127.0.0.1", "localhost", "::1"}:
        return
    raise SynonymDatabaseError("Synonym updates must use HTTPS (localhost HTTP is allowed for testing).")


def validate_synonym_database(payload):
    if not isinstance(payload, dict):
        raise SynonymDatabaseError("The synonym database must be a JSON object.")
    if payload.get("schema_version") != SCHEMA_VERSION:
        raise SynonymDatabaseError("Unsupported synonym database schema version.")

    revision = payload.get("database_revision")
    version = payload.get("dictionary_version")
    dictionaries = payload.get("dictionaries")
    if not isinstance(revision, int) or revision < 1:
        raise SynonymDatabaseError("database_revision must be a positive integer.")
    if not isinstance(version, str) or not version.strip():
        raise SynonymDatabaseError("dictionary_version must be a non-empty string.")
    if not isinstance(dictionaries, dict):
        raise SynonymDatabaseError("dictionaries must be a JSON object.")

    normalized = {
        "schema_version": SCHEMA_VERSION,
        "database_revision": revision,
        "dictionary_version": version.strip(),
        "dictionaries": {},
    }
    for database_name in SUPPORTED_DATABASES:
        mapping = dictionaries.get(database_name)
        if not isinstance(mapping, dict):
            raise SynonymDatabaseError(f"Missing or invalid dictionary: {database_name}.")
        clean_mapping = {}
        alias_owners = {}
        for gene, aliases in mapping.items():
            if not isinstance(gene, str) or not gene.strip():
                raise SynonymDatabaseError(f"Invalid canonical gene in {database_name}.")
            if not isinstance(aliases, list) or not all(
                isinstance(alias, str) and alias.strip() for alias in aliases
            ):
                raise SynonymDatabaseError(
                    f"Aliases for {database_name}/{gene} must be non-empty strings."
                )
            canonical_gene = gene.strip()
            clean_aliases = list(dict.fromkeys(alias.strip() for alias in aliases))
            for alias in [canonical_gene, *clean_aliases]:
                normalized_alias = _normalize_alias(alias)
                previous_gene = alias_owners.get(normalized_alias)
                if previous_gene and previous_gene != canonical_gene:
                    raise SynonymDatabaseError(
                        f"Conflicting alias in {database_name}: {alias!r} belongs to "
                        f"both {previous_gene!r} and {canonical_gene!r}."
                    )
                alias_owners[normalized_alias] = canonical_gene
            clean_mapping[canonical_gene] = clean_aliases
        normalized["dictionaries"][database_name] = clean_mapping
    return normalized


def _load_database_bytes(raw_bytes):
    try:
        payload = json.loads(raw_bytes.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SynonymDatabaseError("The synonym database is not valid UTF-8 JSON.") from exc
    return validate_synonym_database(payload)


def _load_database_path(path):
    return _load_database_bytes(Path(path).read_bytes())


def _database_hash(database):
    canonical = json.dumps(
        database, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def get_bundled_database():
    return _load_database_path(_bundled_database_path())


def get_effective_database():
    bundled = get_bundled_database()
    cached_path = get_cached_database_path()
    if cached_path.is_file():
        try:
            cached = _load_database_path(cached_path)
            if cached["database_revision"] > bundled["database_revision"]:
                return cached, "downloaded"
        except (OSError, SynonymDatabaseError):
            pass
    return bundled, "bundled"


def get_synonym_database_info():
    database, source = get_effective_database()
    return {
        "schema_version": database["schema_version"],
        "database_revision": database["database_revision"],
        "dictionary_version": database["dictionary_version"],
        "source": source,
        "sha256": _database_hash(database),
    }


def get_official_alias_map(org_type):
    database, _ = get_effective_database()
    mapping = database["dictionaries"].get(org_type, {})
    return {gene: aliases.copy() for gene, aliases in mapping.items()}


def _ssl_context():
    try:
        import certifi

        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        return ssl.create_default_context()


def _fetch_bytes(url, max_bytes, timeout=12):
    _validate_remote_url(url)
    request = Request(url, headers={"User-Agent": "Pynnotate-Synonym-Updater/1.0"})
    with urlopen(request, timeout=timeout, context=_ssl_context()) as response:
        content_length = response.headers.get("Content-Length")
        if content_length and int(content_length) > max_bytes:
            raise SynonymDatabaseError("The downloaded synonym file is unexpectedly large.")
        content = response.read(max_bytes + 1)
    if len(content) > max_bytes:
        raise SynonymDatabaseError("The downloaded synonym file is unexpectedly large.")
    return content


def _parse_manifest(raw_bytes, manifest_url):
    try:
        manifest = json.loads(raw_bytes.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SynonymDatabaseError("The synonym update manifest is invalid JSON.") from exc
    if not isinstance(manifest, dict) or manifest.get("schema_version") != SCHEMA_VERSION:
        raise SynonymDatabaseError("Unsupported synonym update manifest.")
    revision = manifest.get("database_revision")
    version = manifest.get("dictionary_version")
    expected_hash = manifest.get("sha256")
    database_url = manifest.get("database_url")
    if not isinstance(revision, int) or revision < 1:
        raise SynonymDatabaseError("Invalid database revision in the update manifest.")
    if not isinstance(version, str) or not version.strip():
        raise SynonymDatabaseError("Invalid dictionary version in the update manifest.")
    if not isinstance(expected_hash, str) or len(expected_hash) != 64:
        raise SynonymDatabaseError("Invalid checksum in the update manifest.")
    if not isinstance(database_url, str) or not database_url.strip():
        raise SynonymDatabaseError("Missing database URL in the update manifest.")
    resolved_url = urljoin(manifest_url, database_url.strip())
    _validate_remote_url(resolved_url)
    return {
        "database_revision": revision,
        "dictionary_version": version.strip(),
        "sha256": expected_hash.lower(),
        "database_url": resolved_url,
    }


def _change_summary(current, candidate):
    genes_added = genes_removed = aliases_added = aliases_removed = 0
    for database_name in SUPPORTED_DATABASES:
        current_map = current["dictionaries"][database_name]
        candidate_map = candidate["dictionaries"][database_name]
        genes_added += len(set(candidate_map) - set(current_map))
        genes_removed += len(set(current_map) - set(candidate_map))
        for gene in set(current_map) | set(candidate_map):
            current_aliases = set(current_map.get(gene, []))
            candidate_aliases = set(candidate_map.get(gene, []))
            aliases_added += len(candidate_aliases - current_aliases)
            aliases_removed += len(current_aliases - candidate_aliases)
    return {
        "genes_added": genes_added,
        "genes_removed": genes_removed,
        "aliases_added": aliases_added,
        "aliases_removed": aliases_removed,
    }


def check_for_synonym_updates(manifest_url=None):
    manifest_url = (
        manifest_url
        or os.environ.get("PYNNOTATE_SYNONYM_MANIFEST_URL", "").strip()
        or DEFAULT_MANIFEST_URL
    )
    _validate_remote_url(manifest_url)
    manifest = _parse_manifest(
        _fetch_bytes(manifest_url, MAX_MANIFEST_BYTES), manifest_url
    )
    current, current_source = get_effective_database()
    if manifest["database_revision"] <= current["database_revision"]:
        return {
            "update_available": False,
            "current_version": current["dictionary_version"],
            "current_revision": current["database_revision"],
            "current_source": current_source,
            "remote_version": manifest["dictionary_version"],
            "remote_revision": manifest["database_revision"],
        }

    database_bytes = _fetch_bytes(manifest["database_url"], MAX_DATABASE_BYTES)
    actual_hash = hashlib.sha256(database_bytes).hexdigest()
    if actual_hash != manifest["sha256"]:
        raise SynonymDatabaseError("The synonym database checksum does not match the manifest.")
    candidate = _load_database_bytes(database_bytes)
    if (
        candidate["database_revision"] != manifest["database_revision"]
        or candidate["dictionary_version"] != manifest["dictionary_version"]
    ):
        raise SynonymDatabaseError("The synonym database version does not match the manifest.")
    return {
        "update_available": True,
        "current_version": current["dictionary_version"],
        "current_revision": current["database_revision"],
        "current_source": current_source,
        "remote_version": candidate["dictionary_version"],
        "remote_revision": candidate["database_revision"],
        "sha256": actual_hash,
        "database_bytes": database_bytes,
        "changes": _change_summary(current, candidate),
    }


def install_synonym_update(update_result):
    if not update_result.get("update_available"):
        raise SynonymDatabaseError("No synonym database update is available.")
    database_bytes = update_result.get("database_bytes")
    if not isinstance(database_bytes, bytes):
        raise SynonymDatabaseError("Missing synonym database update content.")
    actual_hash = hashlib.sha256(database_bytes).hexdigest()
    if actual_hash != update_result.get("sha256"):
        raise SynonymDatabaseError("The synonym database checksum changed before installation.")
    candidate = _load_database_bytes(database_bytes)
    current, _ = get_effective_database()
    if candidate["database_revision"] <= current["database_revision"]:
        raise SynonymDatabaseError("The downloaded synonym database is not newer.")

    data_dir = get_synonym_data_dir()
    data_dir.mkdir(parents=True, exist_ok=True)
    destination = get_cached_database_path()
    file_descriptor, temporary_name = tempfile.mkstemp(
        dir=data_dir, prefix="synonyms_database.", suffix=".tmp"
    )
    try:
        with os.fdopen(file_descriptor, "wb") as handle:
            handle.write(database_bytes)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, destination)
    finally:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
    return get_synonym_database_info()


def restore_bundled_synonym_database():
    cached_path = get_cached_database_path()
    if cached_path.is_file():
        cached_path.unlink()
        return True
    return False


def build_synonym_snapshot(org_type, effective_map, custom_synonyms_applied=False):
    """Build the reproducibility payload saved alongside every result set."""
    info = get_synonym_database_info()
    payload = {
        "schema_version": SCHEMA_VERSION,
        "database_revision": info["database_revision"],
        "dictionary_version": info["dictionary_version"],
        "database_source": info["source"],
        "database_sha256": info["sha256"],
        "genome_type": org_type,
        "custom_synonyms_applied": bool(custom_synonyms_applied),
        "synonyms": {
            gene: list(aliases) for gene, aliases in effective_map.items()
        },
    }
    payload["effective_sha256"] = _database_hash(payload["synonyms"])
    return payload


def export_effective_synonyms(org_type, user_synonyms=None):
    official_map = get_official_alias_map(org_type)
    effective_map = {gene: aliases.copy() for gene, aliases in official_map.items()}
    for gene, aliases in (user_synonyms or {}).items():
        if isinstance(aliases, str):
            aliases = [aliases]
        aliases = list(aliases or [])
        effective_map[gene] = list(dict.fromkeys(effective_map.get(gene, []) + aliases))
    return build_synonym_snapshot(org_type, effective_map, bool(user_synonyms))
