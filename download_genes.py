## modificar sempre o email (linha 22) e a lista de números de acesso (linha 25)

import logging
import traceback
import pandas as pd
from Bio import Entrez, SeqIO
from datetime import datetime

## configurando um arquivo log
date_string = datetime.now().strftime("%d%b%y").lower()
time_string = datetime.now().strftime("%Hh%M").lower()
logging.basicConfig(
    level = logging.INFO,
    format = '%(asctime)s - %(levelname)s - %(message)s',
    handlers = [
        logging.FileHandler(f"download_genes_{date_string}_{time_string}.log"),
        logging.StreamHandler()
    ]
)

## definir o e-mail para o ncbi
Entrez.email = "fernandadesouzacaron@gmail.com"

## lista de números de acesso para fazer download
accession_ids = ["Y16884.3", "MN356140.1", "NC_032362.1", "NC_032363.1", "LC721527.1", "LC715365.1", "LC714913.1", "LC528834.1", "Y18522.2", "AP009194.1", "AP009193.1", "AP009191.1", "AP009189.1", "LC541480.1", "LC541479.1", "LC541477.1", "LC541476.1", "LC541475.1", "LC541474.1", "LC541473.1", "LC541472.1", "LC541470.1", "LC541467.1", "LC541462.1", "LC541461.1", "LC541460.1", "LC541458.1", "LC541457.1", "LC541456.1", "LC541455.1", "LC541453.1", "LC541452.1", "LC541449.1", "LC541448.1", "LC541447.1", "LC541445.1", "LC541443.1", "LC541438.1", "LC541437.1", "LC541432.1", "LC541431.1", "LC541430.1", "LC541429.1", "LC541469.1", "LC541463.1", "LC541450.1", "LC541446.1", "LC541444.1", "LC541441.1", "LC541440.1", "LC532223.1", "AB164627.1", "AB164626.1", "AB164625.1", "AB164624.1", "AB164622.1", "AM889141.1", "AM889140.1", "AM889139.1", "AB073301.1", "AP005595.1", "AP010821.1", "AP003325.1", "AP006741.1", "AP003324.1", "KM078794.1", "AP010797.1", "AP009043.1", "AP009042.1", "AB026818.1", "AB026193.1", "AP010823.1", "AP010822.1", "AM237310.1", "AP008238.1", "AP008239.1", "NC_023777.1", "NC_000877.1", "MH041488.1", "MH041485.1", "MH041486.1", "NC_018548.1", "NC_000879.1", "MT341548.1", "NC_082066.1", "NC_039536.1", "NC_018033.1", "NC_000880.1", "OP718272.1", "KX245139.1", "KX245135.1", "EU755254.1", "NC_086745.1", "MT319112.1", "NC_071761.1", "NC_071238.1", "NC_071237.1", "NC_071234.1", "NC_059907.1", "NC_056091.1", "NC_051886.1", "NC_050289.1", "NC_046877.1", "NC_045376.1", "NC_045370.1", "NC_045366.1", "NC_039538.1", "NC_036760.1", "NC_030604.1", "NC_030603.1", "NC_029359.1", "NC_029328.1", "NC_028414.1", "NC_028082.1", "NC_027657.1", "NC_027241.1", "NC_027096.1", "NC_027095.1", "NC_026911.1", "NC_026715.1", "NC_025742.1", "NC_024553.1", "NC_023264.1", "NC_021105.1", "NC_020425.1", "NC_015482.1", "NC_015234.1", "NC_015087.1", "NC_015085.1", "NC_014455.1", "NC_014341.1", "NC_014046.1", "NC_003128.3", "OP936976.1", "OP936974.1", "OP936970.1", "OP936964.1", "OL639163.1", "MW364567.1", "MZ048030.1", "MT366880.1", "MN067867.1", "MK377247.1", "MF071218.1", "KU886092.1", "KC836121.1", "MF431745.1", "EU165707.1", "EU165706.1", "FJ799728.1", "KU221053.1", "KP137624.1", "JQ282801.1", "KU237286.1", "PP898113.1", "PP576342.1", "OP998296.1", "OQ561760.1", "PP357441.1", "NC_086816.1", "PP174885.1", "NC_084285.1", "NC_083117.1", "OR120383.1", "MN356405.1", "NC_081048.1", "NC_077264.1", "NC_061952.1", "NC_059847.1", "NC_064933.1", "NC_059854.1", "NC_059853.1", "NC_059848.1", "NC_059801.1", "NC_066800.1", "NC_062945.1", "NC_062944.1", "NC_072120.1", "NC_070231.1", "NC_069920.1", "NC_062819.1", "NC_060992.1", "NC_058320.1", "NC_057528.1", "NC_057527.1", "NC_057295.1", "NC_056389.1", "NC_056273.1", "NC_056257.1", "NC_054361.1", "NC_054275.1", "NC_054153.1", "NC_053917.1", "NC_053867.1", "NC_051558.1", "NC_051538.1", "NC_051513.1", "NC_050973.1", "NC_050680.1", "NC_050393.1", "NC_050299.1", "NC_050293.1", "NC_050291.1", "NC_050290.1", "NC_050287.1", "NC_050286.1", "NC_050278.1", "NC_050259.1", "NC_050049.1", "NC_046765.1", "NC_046029.1", "NC_045517.1", "NC_045401.1", "NC_045381.1", "NC_045377.1", "NC_045375.1", "NC_045374.1", "NC_045373.1", "NC_045372.1", "NC_045371.1", "NC_045369.1", "NC_045368.1", "NC_045367.1", "NC_045364.1", "NC_045282.1", "NC_045181.1", "NC_045042.1", "NC_044743.1", "NC_044648.1", "NC_043899.1", "NC_042414.1", "NC_042191.1", "NC_041422.1", "NC_041257.1", "NC_041197.1", "NC_041141.1", "NC_041121.1", "NC_041111.1", "NC_041109.1", "NC_041095.1", "NC_040991.1", "NC_040990.1", "NC_040956.1", "NC_040850.1", "NC_040290.1", "NC_040142.1", "NC_039950.1", "NC_039888.1", "NC_039842.1", "NC_039770.1", "NC_039736.1", "NC_039396.1", "NC_038201.1", "NC_038195.1", "NC_037692.1", "NC_037486.1", "NC_037464.1", "NC_036016.1", "NC_035877.1", "NC_035806.1", "NC_035801.1", "NC_035868.1", "NC_035747.1", "NC_035746.1", "NC_035673.1", "NC_035658.1", "NC_035568.1", "NC_034679.1", "NC_034642.1", "NC_034373.1", "NC_034353.1", "NC_034296.1", "NC_034278.1", "NC_034237.1", "NC_033974.1", "NC_033966.1", "NC_033973.1", "NC_033967.1", "NC_033536.1", "NC_033414.1", "NC_033413.1", "NC_033406.1", "NC_033336.1", "NC_032059.1", "NC_032058.1", "NC_031869.1", "NC_031866.1", "NC_031830.1", "NC_031819.1", "NC_031867.1", "NC_031865.1", "NC_031863.1", "NC_031447.1", "NC_031353.1", "NC_031350.1", "NC_030771.1", "NC_030588.1", "NC_030585.1", "NC_030507.1", "NC_030288.1", "NC_030287.1", "NC_030286.1", "NC_030285.1", "NC_029862.1", "NC_029769.1", "NC_029482.1", "NC_029475.1", "NC_029462.1", "NC_029402.1", "NC_029377.1", "NC_029344.1", "NC_029322.1", "NC_029321.1", "NC_029319.1", "NC_029188.1", "NC_029161.1", "NC_029189.1", "NC_028545.1", "NC_028437.1", "NC_028436.1", "NC_028408.1", "NC_028404.1", "NC_028346.1", "NC_028333.1", "NC_028195.1", "NC_028194.1", "NC_028193.1", "NC_028187.1", "NC_028186.1", "NC_028179.1", "NC_028039.1", "NC_028038.1", "NC_028037.1", "NC_028036.1", "NC_028020.1", "NC_028019.1", "NC_027942.1", "NC_027936.1", "NC_027934.1", "NC_027933.1", "NC_027606.1", "NC_027504.1", "NC_027496.1", "NC_027420.1", "NC_027285.1", "NC_027260.1", "NC_027284.1", "NC_027279.1", "NC_026793.1", "NC_026701.1", "NC_026461.1", "NC_026459.1", "NC_026223.1", "NC_026082.1", "NC_026068.1", "NC_026067.1", "NC_026031.1", "NC_026029.1", "NC_025925.1", "NC_025924.1", "NC_025923.1", "NC_025922.1", "NC_025921.1", "NC_025920.1", "NC_025919.1", "NC_025917.1", "NC_025916.1", "NC_025900.1", "NC_025786.1", "NC_025654.1", "NC_025649.1", "NC_025637.1", "NC_025579.1", "NC_025556.1", "NC_025521.1", "NC_025318.1", "NC_024922.1", "NC_024750.1", "NC_024730.1", "NC_024726.1", "NC_024640.1", "NC_024620.1", "NC_024619.1", "NC_024618.1", "NC_024617.1", "NC_024616.1", "NC_024615.1", "NC_024602.1", "NC_024595.1", "NC_024554.1", "NC_024552.1", "NC_024539.1", "NC_024533.1", "NC_024267.1", "NC_024266.1", "NC_024257.1", "NC_024156.1", "NC_024107.1", "NC_024069.1", "NC_024048.1", "NC_024068.1", "NC_023982.1", "NC_022817.1", "NC_022418.1", "NC_021771.1", "NC_021641.1", "NC_021368.1", "NC_020605.1", "NC_020604.1", "NC_020603.1", "NC_020602.1", "NC_020601.1", "NC_020600.1", "NC_020599.1", "NC_020598.1", "NC_020597.1", "NC_020596.1", "NC_020582.1", "NC_020581.1", "NC_020580.1", "NC_020579.1", "NC_020577.1", "NC_020576.1", "NC_020574.1", "NC_020570.1", "NC_020569.1", "NC_020427.1", "NC_020426.1", "NC_019804.1", "NC_018034.1", "NC_016723.1", "NC_016679.1", "NC_015898.1", "NC_015897.1", "NC_015887.1", "NC_015810.1", "NC_015802.1", "NC_015613.1", "NC_015526.1", "NC_015237.1", "NC_015236.1", "NC_015233.1", "NC_015232.1", "NC_015198.1", "NC_015114.1", "NC_015074.1", "NC_014879.1", "NC_013979.1", "NC_013619.1", "NC_012844.1", "NC_012843.1", "NC_012453.1", "NC_009134.1", "NC_007172.2", "NC_002778.2", "NC_002772.2", "NC_002781.3", "NC_002782.2", "NC_002783.2", "NC_003713.2", "NC_003712.2", "NC_002785.1", "NC_002784.1", "EU376027.1", "ON920867.1", "MW417350.1", "OK542103.1", "MN248536.1", "MT341558.1", "OK662584.1", "MW489467.1", "MZ677204.1", "MZ681908.1", "MW880933.1", "MW574395.1", "MW574394.1", "MW574393.1", "MW574391.1", "MW574390.1", "MW574389.1", "MW574387.1", "MW574386.1", "MW574385.1", "MW574384.1", "MW574382.1", "MW574381.1", "MW574380.1", "MW574379.1", "MW574377.1", "MW574376.1", "MW574375.1", "MW574374.1", "MW574373.1", "MW574372.1", "MW574370.1", "MW574368.1", "MW574367.1", "MW574365.1", "MW574363.1", "MW574362.1", "MW574361.1", "MW574360.1", "MW574359.1", "MW574358.1", "MW574357.1", "MW574356.1", "MW574355.1", "MW574354.1", "MW574353.1", "MW574352.1", "MW574351.1", "MW574349.1", "MZ337397.1", "MT792356.1", "MW861756.1", "MW338659.1", "MG833030.1", "MW495245.1", "MN206975.1", "MT427586.1", "MN833781.1", "MT471263.1", "MF683387.1", "MT181966.1", "MN857545.1", "MN481404.1", "MN217252.1", "MK992912.1", "MN515396.1", "MK051002.1", "MK408609.1", "MN047457.1", "MK940810.1", "MK714020.1", "KX977449.1", "MN125374.1", "MK905885.1", "MK342599.1", "MN122908.1", "MH433599.1", "MG583885.1", "MH133972.1", "MH133971.1", "MH133970.1", "MH133969.1", "MH133968.1", "EU417810.1", "MF977813.1", "MG596878.1", "KY349099.1", "MF784450.1", "KY767670.1", "KY888681.1", "KX925978.1", "KY419885.1", "KF946546.1", "KX809695.1", "KX902239.1", "EU417812.1", "EU417811.1", "HQ890328.1", "HQ915865.1", "KT427463.1", "KU884610.1", "KP684122.1", "KP019940.1", "KP995437.1", "KT633399.1", "KR072661.1", "KJ680300.1", "KP313718.1", "KF951092.1", "KF951088.1", "KM873665.1", "JX524614.1", "JX215256.1", "KJ190959.1", "KJ190954.1", "KJ190950.1", "KJ190949.1", "KJ190957.1", "KF437906.1", "KC953095.1", "KC936100.1", "KJ834096.1", "KF293721.1", "KJ147475.1", "KF682364.1", "JX524615.1", "JQ782214.1", "JQ003192.1", "FJ769845.1", "OR731193.2", "MT547766.1", "NC_087758.1", "OR750915.1", "NC_083486.1", "NC_080969.1", "NC_080911.1", "NC_080910.1", "OR030351.1", "OR030350.1", "OP609734.1", "NC_072569.1", "NC_013806.1", "NC_068694.1", "NC_068693.1", "NC_068692.1", "NC_068688.1", "NC_068687.1", "NC_060467.1", "NC_057294.1", "NC_052013.1", "NC_052011.1", "NC_051552.1", "NC_050688.1", "NC_045378.1", "NC_041666.1", "NC_041665.1", "NC_041663.1", "NC_041662.1", "NC_041661.1", "NC_041660.1", "NC_041659.1", "NC_041165.1", "NC_041118.1", "NC_039818.1", "NC_039537.1", "NC_038220.1", "NC_038152.1", "NC_037702.1", "NC_037154.1", "NC_037153.1", "NC_036337.1", "NC_036297.1", "NC_036050.1", "NC_034933.1", "NC_034316.1", "NC_033338.1", "NC_032725.1", "NC_031897.1", "NC_031845.1", "NC_029837.1", "NC_029384.1", "NC_029383.1", "NC_029360.1", "NC_029340.1", "NC_028441.1", "NC_028178.1", "NC_028177.1", "NC_028163.1", "NC_027848.1", "NC_027847.1", "NC_027846.1", "NC_027845.1", "NC_027844.1", "NC_027843.1", "NC_027842.1", "NC_027841.1", "NC_027840.1", "NC_027839.1", "NC_027251.1", "NC_027231.1", "NC_027230.1", "NC_026548.1", "NC_026547.1", "NC_025927.1", "NC_025926.1", "NC_024925.1", "NC_024924.1", "NC_024682.1", "NC_024673.1", "NC_023940.1", "NC_023939.1", "NC_022840.1", "NC_022839.1", "NC_022452.1", "NC_022150.1", "NC_021445.1", "NC_021408.1", "NC_019668.1", "NC_019667.1", "NC_019664.1", "NC_010195.2", "NC_013483.2", "NC_011307.1", "NC_010089.1", "NC_010095.1", "NC_010094.1", "NC_010091.1", "NC_008549.1", "NC_008551.1", "NC_008547.1", "NC_008548.1", "NC_008550.1", "NC_008546.1", "NC_007691.1", "NC_007174.1", "NC_007007.1", "NC_007011.1", "NC_005932.1", "NC_005933.1", "NC_004538.1", "ON986363.1", "OK539292.1", "MZ128784.1", "OM063155.2", "MN431594.1", "MW623728.1", "MW623686.1", "MG883743.1", "MN125373.1", "MK820678.1", "MH700653.1", "MH700651.1", "MH700645.1", "MH700633.1", "MH700631.1", "MG681082.1", "LC099103.1", "KU361806.1", "KU356676.1", "KU356673.1", "KY994613.1", "KY994608.1", "KY994606.1", "KY994604.1", "KY994598.1", "KY994597.1", "KY994595.1", "KY994594.1", "KY994593.1", "KY994592.1", "KY994589.1", "KY994588.1", "KY994587.1", "KY994585.1", "KY994581.1", "KT340629.1", "KX592585.1", "KY419888.1", "KT901459.1", "DQ453514.1", "KM251414.1", "KJ631624.1", "KT934323.1", "KM374660.1", "KM374659.1", "KM374655.1", "KJ716444.1", "NC_053923.1", "KX245144.1", "KU094576.1", "MN065674.1", "MK714019.1", "MH041272.1", "KP893332.1", "KX902241.1", "KX902237.1", "KX189345.1", "AY016012.1", "AY016011.1", "AY016010.1", "KJ680303.1", "KJ914547.1", "BK062963.1", "NC_059932.1"]

## dicionário de sinônimos
synonym_dict = {
    "ATP6": ["ATP6", "atp6", "ATPase6", "ATPase 6", "ATP synthase 6", 
             "ATP synthase subunit 6", "ATP synthase F0 subunit 6", 
             "ATPase subunit 6", "MT-ATP6", "mt-Atp6", "mt-atp6",
             "F0-ATP synthase subunit6", "atpase6", "atpase 6", "Atp6",
             "ATP sythase subunit 6", "AT6", "MTATP6"],
    "ATP8": ["ATP8", "atp8", "ATPase8", "ATPase 8", "ATP synthase 8", 
             "ATP synthase subunit 8", "ATP synthase F0 subunit 8", 
             "ATPase subunit 8", "MT-ATP8", "mt-Atp8", "mt-atp8", "Atp8",
             "F0-ATP synthase subunit8", "atpase8", "atpase 8",
             "ATP sythase subunit 8", "AT8", "MTATP8"],
    "COI": ["cytochrome c oxidase subunit 1", "COI", "COX1", "cox1", "CO1", 
            "COXI", "cytochrome c oxidase subunit I", "COX-I", "coi", "MT-CO1",
            "mt-Co1", "mt-co1", "cytochrome oxidase c subunit 1",
            "cytochrome oxidase subunit 1", "cytochrome oxidase subunit1",
            "Cytochrome c oxidase subunit1", "CO I", "coi", "co1", "coI",
            "coxI", "cytochrome oxidase I", "cytochrome oxidase subunit I",
            "cox I", "Cox1"],
    "COII": ["cytochrome c oxidase subunit 2", "COII", "COX2", "cox2", "COXII",
             "CO2", "cytochrome c oxidase subunit II", "COX-II", "MT-CO2", 
             "mt-Co2", "mt-co2", "cytochrome oxidase subunit 2", "coxII",
             "cytochrome oxidase subunit 2", "cytochrome oxidase subunit2",
             "Cytochrome c oxidase subunit2", "CO II", "coii", "co2", "coII",
             "cytochrome oxidase II", "cytochrome oxidase subunit II", 
             "cox II", "Cox2"],
    "COIII": ["cytochrome c oxidase subunit 3", "COIII", "COX3", "cox3", "CO3", 
              "COXIII", "cytochrome c oxidase subunit III", "COX-III", "MT-CO3",
               "mt-Co3", "mt-co3", "cytochrome oxidase c subunit 3", "coxIII",
               "cytochrome oxidase subunit 3", "cytochrome oxidase subunit3",
               "Cytochrome c oxidase subunit3", "CO III", "coiii",
               "cytochrome c oxidase subunit3", "co3", "coIII",
               "cytochrome oxidase III", "cytochrome oxidase subunit III", 
               "CO3 subunit 3", "cox III", "Cox3"],
    "CYTB": ["CYTB", "cytb", "cob", "Cyt B", "Cyt b", "Cytb", "cyt b", "Cb", 
             "cytochrome b", "CYB", "cytB", "cyb", "MT-CYB", "mt-cyb", "cyt-B",
             "mt-Cytb", "cytochorome b", "Cytochrome b", "ctyb", "COB"],
    "ND1": ["ND1", "NADH1", "NADH dehydrogenase subunit 1", "nad1", "nd1", 
            "MT-ND1", "mt-Nd1", "mt-nd1", "NADH-1", "MTND1", "nadh1",
            "NADH dehydrogenase subunit1", "NAD1", "NADH subunit 1", 
            "NADH dehydrogenase 1", "NADH dehydrogenase subunit I"],
    "ND2": ["ND2", "NADH2", "NADH dehydrogenase subunit 2", "nad2", "nd2", 
            "NADH dehydrogenase subunit II", "NADH subunit 2", "MT-ND2", 
            "mt-Nd2", "mt-nd2", "NAD2", "NADH dehydrogenase subunit2", "nadh2",
            "NAD2", "NADH subunit 2"],
    "ND3": ["ND3", "NADH3", "NADH dehydrogenase subunit 3", "nad3", "nd3", 
            "MT-ND3", "mt-Nd3", "mt-nd3", "NADH dehydrogenase subunit3", 
            "nadh3", "NAD3", "NADH subunit 3", 
            "NADH dehydrogenase subunit III"],
    "ND4": ["ND4", "NADH4", "NADH dehydrogenase subunit 4", "nad4", "nd4", 
            "MT-ND4", "mt-Nd4", "mt-nd4", "NADH dehydrogenase subunit4", 
            "nadh4", "4", "NAD4", "NADH subunit 4",
            "NADH dehydrogenase subunit IV"],
    "ND4L": ["ND4L", "NADH4L", "NADH dehydrogenase subunit 4L", "nad4l", "nd4l",
             "nd4L", "MT-ND4L", "mt-Nd4l", "mt-nd4l", "nad4L", "nadh4L",
             "NADH dehydrogenase subunit 4 L", "NADH dehydrogenase subunit4L",
             "NAD4L", "NADH subunit 4L", "ND4l",
             "NADH dehydrogenase subunit IV L"],
    "ND5": ["ND5", "NADH5", "NADH dehydrogenase subunit 5", "nad5", "nd5", 
            "MT-ND5", "mt-Nd5", "mt-nd5", "nadh5", "nadh5", "NAD5",
            "NADH dehydrogenase subunit5", "NADH subunit 5",
            "NADH dehydrogenase subunit V"],
    "ND6": ["ND6", "NADH6", "NADH dehydrogenase subunit 6", "nad6", "nd6", 
            "MT-ND6", "mt-Nd6", "mt-nd6", "NADH dehydrogenase subunit6",
            "nadh6", "NAD6", "NADH subunit 6",
            "NADH dehydrogenase subunit VI"],
    "12S": ["12S ribosomal RNA", "s-rRNA", "small subunit ribosomal RNA", "12S",
            "12S rrn", "12SrRNA", "MTRNR1", "mt-Rnr1", "mt-rnr1", "mtrnr1", 
            "MT-RNR1", "SSU", "ssu", "rrn12", "ssu rRNA", "12 S ribosomal RNA",
            "small ribosomal RNA subunit RNA", "12S rRNA", "rRNA-12S",
            "12S ribosormal RNA", "12S-rRNA"],
    "16S": ["16S ribosomal RNA", "l-rRNA", "large subunit ribosomal RNA", "16S",
            "16rrn", "16S rrn", "16S rRNA", "16Srrn", "16SrRNA", "lsu", "LSU", 
            "lsu rRNA", "MTRNR2", "mt-Rnr2", "mt-rnr2", "MT-RNR2", "rrn16",
            "large ribosomal RNA subunit RNA", "16 S ribosomal RNA", "l-RNA",
            "16S-rRNA"],
    "tRNA_Ala": ["tRNA-Ala", "trnA", "trnA-ugc", "trnA TGC"],
    "tRNA_Arg": ["tRNA-Arg", "trnR", "trnR-ucg", "trnR TCG"],
    "tRNA_Asn": ["tRNA-Asn", "trnN", "trnN-guu", "trnN GTT"],
    "tRNA_Asp": ["tRNA-Asp", "trnD", "trnD-guc", "trnD GTC"],
    "tRNA_Cys": ["tRNA-Cys", "trnC", "trnC-gca", "trnC GCA"],
    "tRNA_Gln": ["tRNA-Gln", "trnQ", "trnQ-uug", "trnQ TTG"],
    "tRNA_Glu": ["tRNA-Glu", "trnE", "trnE-uuc", "trnE TTC"],
    "tRNA_Gly": ["tRNA-Gly", "trnG", "trnG-ucc", "trnG TCC"],
    "tRNA_His": ["tRNA-His", "trnH", "trnH-gug", "trnH GTG"],
    "tRNA_Ile": ["tRNA-Ile", "trnI", "trnI-gau", "trnI GAT"],
    "tRNA_Leu": ["tRNA-Leu", "trnL", "trnL-uag", "trnL TAG", "tRNA-Leu (CUN)",
                 "tRNA-Leu (UUR)", "tRNA-Leu(CUN)", "tRNA-Leu(UUR)"],
    "tRNA_Lys": ["tRNA-Lys", "trnK", "trnK-uuu", "trnK TTT"],
    "tRNA_Met": ["tRNA-Met", "trnM", "trnM-cau", "trnM CAT"],
    "tRNA_Phe": ["tRNA-Phe", "trnF-gaa", "trnF GAA"],
    "tRNA_Pro": ["tRNA-Pro", "trnP-ugg", "trnP TGG", "proline tRNA"],
    "tRNA_Ser": ["tRNA-Ser", "trnS", "trnS-uga", "trnS GCT", "tRNA-Ser (UCN)",
                 "tRNA-Ser (AGY)", "tRNA-Ser(UCN)", "tRNA-Ser(AGY)"],
    "tRNA_Thr": ["tRNA-Thr", "trnT-ugu", "trnT TGT", "threonine tRNA"],
    "tRNA_Trp": ["tRNA-Trp", "trnW", "trnW-uca", "trnW TCA"],
    "tRNA_Tyr": ["tRNA-Tyr", "trnY", "trnY-gua", "trnY GTA"],
    "tRNA_Val": ["tRNA-Val", "trnV", "trnV-uac", "trnV TAC"]
}

## função para encontrar o nome geral com base no sinônimo
def get_official_name(feature_name):
    for official_name, synonyms in synonym_dict.items():
        if feature_name in synonyms:
            return official_name
    return None

## função para baixar e salvar as sequências das features
def download_and_save_genes(genbank_id):
    
    ## baixar o arquivo do GenBank
    handle = Entrez.efetch(db = "nucleotide", id = genbank_id, rettype = "gb",
                           retmode = "text")
    ## lendo o arquivo baixado no formato genbank
    record = SeqIO.read(handle, "genbank")
    ## fecha conexão
    handle.close()

    ## pegar o nome da espécie e substitui espaços "_"
    species_name = record.annotations.get("organism", "unknown_species").replace(" ", "_")

    ## printando qual espécies está sendo feito o download
    logging.info(f"Obtendo sequências de {species_name} ({genbank_id})...")

    # criar um objeto para armazenar as features encontradas
    features_found = {}

    ## criar um objeto para armazenar as features NÃO encontradas
    missing_features = set(synonym_dict.keys())

    ## dicionário para armazenar as localizações dos tRNAs com mais de uma cópia
    tRNA_locations = {
        "tRNA_Leu": [],
        "tRNA_Ser": []
    }

    ## separar as features anotadas (CDS, rRNA, tRNA, D-loop, control_region)
    for feature in record.features:
        
        ## pular se a feature for de um desses tipos
        if feature.type in ["gene", "source"]:  
            continue

        ## se CDS, rRNA, tRNA pegar o nome da feature do campo "gene" ou
        ## "product"
        if feature.type == "CDS":
            gene_name = feature.qualifiers.get("gene", feature.qualifiers.get("product", ["unknown_gene"]))[0]
            feature_type = "CDS"
        elif feature.type == "rRNA":
            gene_name = feature.qualifiers.get("product", feature.qualifiers.get("gene", ["unknown_rRNA"]))[0]
            feature_type = "rRNA"
        elif feature.type == "tRNA":
            gene_name = feature.qualifiers.get("product", feature.qualifiers.get("gene", ["unknown_tRNA"]))[0]
            feature_type = "tRNA"
        else:
            continue

        ## verificar se o nome da feature é um sinônimo conhecido
        official_name = get_official_name(gene_name)

        ## se o nome não for reconhecido, pular a feature
        if not official_name:
            ## printando que não foi encontrado um sinônimo no nosso 
            ## dicionárias para a feature em questão
            logging.error(f"SINÔNIMO desconhecido para {gene_name}, pulando.")
            continue

        ## pegar a sequência da feature
        feature_seq = feature.extract(record.seq)

        ## criar nome para essa sequência (espécies_accession)
        sequence_id = f">{species_name}_{genbank_id}"
        
        ## nome do arquivo baseado no nome oficial da feature com base no
        ## nosso dicionário

        ## primeiro pegar data para nome do arquivo
        date_string = datetime.now().strftime("%d%b%y").lower()

        ## fazer condição para armazenamento do arquivo para considerar 
        ## arquivos separados para tRNAs com mais de uma cópia
        
        ## se for um tRNA e tiver mais de uma cópia, guarda a localização da
        ## feature e a sequência em um objeto separado
        if official_name in tRNA_locations:
            ## pegar localização da feature
            start = int(feature.location.start)
            ## guardar sequência da feature com sua localização naquele 
            ## dicionário que fizemos anteriormente
            tRNA_locations[official_name].append((start, feature_seq))
        ## senão, apenas guarda a feature no arquivo correspondente ao 
        ## gene/rRNA/tRNA em questão
        else:
            ## para outras features, usar um único arquivo
            filename = f"{official_name}_{date_string}.fasta"
            
            ## abrir o arquivo no modo de adição e escrever a sequência
            with open(filename, "a") as output_handle:
                output_handle.write(f"{sequence_id}\n{feature_seq}\n")

        ## remover a feature encontrada da lista de features ausentes
        missing_features.discard(official_name)

        ## armazenar o número de acesso da feature encontrada no dicionário
        features_found[official_name] = genbank_id

    ## depois de analisar todas as features, temos que processar as 
    ## localizações dos tRNAs com mais de uma cópia encontrados
    for tRNA, locations in tRNA_locations.items():
        if len(locations) > 0:
            ## se houver só uma cópia do tRNA, salvar no primeiro arquivo 
            ## do tRNA mesmo
            if len(locations) == 1:
                start, feature_seq = locations[0]
                filename = f"{tRNA}_part1_{date_string}.fasta"
                with open(filename, "a") as output_handle:
                    output_handle.write(f"{sequence_id}\n{feature_seq}\n")
            ## se houver mais de uma cópia, salvar o que aparece primeiro
            ## no primeiro arquivo e o que aparece por segundo, no segundo
            ## arquivo (não é ideal...)
            elif len(locations) == 2:
                start1, seq1 = locations[0]
                start2, seq2 = locations[1]

                if start1 < start2:
                    filename1 = f"{tRNA}_part1_{date_string}.fasta"
                    with open(filename1, "a") as output_handle:
                        output_handle.write(f"{sequence_id}\n{seq1}\n")
                    
                    filename2 = f"{tRNA}_part2_{date_string}.fasta"
                    with open(filename2, "a") as output_handle:
                        output_handle.write(f"{sequence_id}\n{seq2}\n")
                else:
                    filename1 = f"{tRNA}_part1_{date_string}.fasta"
                    with open(filename1, "a") as output_handle:
                        output_handle.write(f"{sequence_id}\n{seq2}\n")
                    
                    filename2 = f"{tRNA}_part2_{date_string}.fasta"
                    with open(filename2, "a") as output_handle:
                        output_handle.write(f"{sequence_id}\n{seq1}\n")

    ## se sobrou alguma feature não encontrada para a espécies, emitir aviso
    for missing_feature in missing_features:
        logging.warning(f"Aviso: '{missing_feature}' não encontrada para {species_name} ({genbank_id})")

    return species_name, features_found
    ## fim do loop das features

## criando tabela para armazenar os números de acesso
data_accession = []

## criado um conjunto para rastrear os nomes das features para colocar na tabela
feature_names = set(synonym_dict.keys())

## loop para processar todos os números de acesso
for accession in accession_ids:
    ## usar o nome da espécie como chave no dicionário
    species_name = None  
    try:
        species_name, features_found = download_and_save_genes(accession)
        ## linha na tabela com o nome da espécie e o accession
        row_data = [species_name, accession]

        ## preencher as colunas com os dados das features
        for official_name in feature_names:
            ## se a feature foi encontrada, adicionar o número de acesso
            if official_name in features_found:
                row_data.append(features_found[official_name])
            else:
                row_data.append("NA")  ## senão, adicionar NA

        ## adicionar a linha à lista de dados da tabela
        data_accession.append(row_data)

    except Exception as e:
        error_info = traceback.format_exc()
        logging.error(f"ERRO ao processar {accession}: {e}")
        logging.error("Detalhes do erro:", error_info)

## criando um dataframe usando o dicionário data_accession
feature_names_prefix = [f"Acc_{name}" for name in feature_names]
data_accession_df = pd.DataFrame(data_accession, 
                                 columns = ["Sequence_names", 
                                 "Mitogenome_accession"] + 
                                 feature_names_prefix)

## salvando a tabela em um arquivo csv
data_accession_df.to_csv(f"spp_accession_{date_string}.csv",
                         index = False)


