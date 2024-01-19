// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - isis agora lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>

//! This module contains backend-specific constant values, such as the 64-bit limbs of curve constants.

use super::field::FieldElement51;
use super::scalar::Scalar52;
use crate::edwards::EdwardsPoint;

#[cfg(feature = "precomputed-tables")]
use crate::{
    backend::serial::curve_models::AffineNielsPoint,
    window::NafLookupTable8,
};

/// Edwards `d` value, equal to `-121665/121666 mod p`.
pub(crate) const EDWARDS_D: FieldElement51 = FieldElement51::from_limbs([
    929955233495203,
    466365720129213,
    1662059464998953,
    2033849074728123,
    1442794654840575,
]);

/// Edwards `2*d` value, equal to `2*(-121665/121666) mod p`.
pub(crate) const EDWARDS_D2: FieldElement51 = FieldElement51::from_limbs([
    1859910466990425,
    932731440258426,
    1072319116312658,
    1815898335770999,
    633789495995903,
]);

/// Precomputed value of one of the square roots of -1 (mod p)
pub(crate) const SQRT_M1: FieldElement51 = FieldElement51::from_limbs([
    1718705420411056,
    234908883556509,
    2233514472574048,
    2117202627021982,
    765476049583133,
]);

/// `APLUS2_OVER_FOUR` is (A+2)/4. (This is used internally within the Montgomery ladder.)
pub(crate) const APLUS2_OVER_FOUR: FieldElement51 =
    FieldElement51::from_limbs([121666, 0, 0, 0, 0]);

/// `MONTGOMERY_A` is equal to 486662, which is a constant of the curve equation
/// for Curve25519 in its Montgomery form. (This is used internally within the
/// Elligator map.)
pub(crate) const MONTGOMERY_A: FieldElement51 = FieldElement51::from_limbs([486662, 0, 0, 0, 0]);

/// `MONTGOMERY_A_NEG` is equal to -486662. (This is used internally within the
/// Elligator map.)
pub(crate) const MONTGOMERY_A_NEG: FieldElement51 = FieldElement51::from_limbs([
    2251799813198567,
    2251799813685247,
    2251799813685247,
    2251799813685247,
    2251799813685247,
]);

/// `L` is the order of base point, i.e. 2^252 + 27742317777372353535851937790883648493
pub(crate) const L: Scalar52 = Scalar52([
    0x0002631a5cf5d3ed,
    0x000dea2f79cd6581,
    0x000000000014def9,
    0x0000000000000000,
    0x0000100000000000,
]);

/// `L` * `LFACTOR` = -1 (mod 2^52)
pub(crate) const LFACTOR: u64 = 0x51da312547e1b;

/// `R` = R % L where R = 2^260
pub(crate) const R: Scalar52 = Scalar52([
    0x000f48bd6721e6ed,
    0x0003bab5ac67e45a,
    0x000fffffeb35e51b,
    0x000fffffffffffff,
    0x00000fffffffffff,
]);

/// `RR` = (R^2) % L where R = 2^260
pub(crate) const RR: Scalar52 = Scalar52([
    0x0009d265e952d13b,
    0x000d63c715bea69f,
    0x0005be65cb687604,
    0x0003dceec73d217f,
    0x000009411b7c309a,
]);


/// The 8-torsion subgroup \\(\mathcal E \[8\]\\).
///
/// In the case of Curve25519, it is cyclic; the \\(i\\)-th element of
/// the array is \\(\[i\]P\\), where \\(P\\) is a point of order \\(8\\)
/// generating \\(\mathcal E\[8\]\\).
///
/// Thus \\(\mathcal E\[4\]\\) is the points indexed by `0,2,4,6`, and
/// \\(\mathcal E\[2\]\\) is the points indexed by `0,4`.
pub(crate) const EIGHT_TORSION: [EdwardsPoint; 8] = EIGHT_TORSION_INNER_DOC_HIDDEN;

/// Inner item used to hide limb constants from cargo doc output.
#[doc(hidden)]
pub(crate) const EIGHT_TORSION_INNER_DOC_HIDDEN: [EdwardsPoint; 8] = [
    EdwardsPoint {
        X: FieldElement51::from_limbs([0, 0, 0, 0, 0]),
        Y: FieldElement51::from_limbs([1, 0, 0, 0, 0]),
        Z: FieldElement51::from_limbs([1, 0, 0, 0, 0]),
        T: FieldElement51::from_limbs([0, 0, 0, 0, 0]),
    },
    EdwardsPoint {
        X: FieldElement51::from_limbs([
            358744748052810,
            1691584618240980,
            977650209285361,
            1429865912637724,
            560044844278676,
        ]),
        Y: FieldElement51::from_limbs([
            84926274344903,
            473620666599931,
            365590438845504,
            1028470286882429,
            2146499180330972,
        ]),
        Z: FieldElement51::from_limbs([1, 0, 0, 0, 0]),
        T: FieldElement51::from_limbs([
            1448326834587521,
            1857896831960481,
            1093722731865333,
            1677408490711241,
            1915505153018406,
        ]),
    },
    EdwardsPoint {
        X: FieldElement51::from_limbs([
            533094393274173,
            2016890930128738,
            18285341111199,
            134597186663265,
            1486323764102114,
        ]),
        Y: FieldElement51::from_limbs([0, 0, 0, 0, 0]),
        Z: FieldElement51::from_limbs([1, 0, 0, 0, 0]),
        T: FieldElement51::from_limbs([0, 0, 0, 0, 0]),
    },
    EdwardsPoint {
        X: FieldElement51::from_limbs([
            358744748052810,
            1691584618240980,
            977650209285361,
            1429865912637724,
            560044844278676,
        ]),
        Y: FieldElement51::from_limbs([
            2166873539340326,
            1778179147085316,
            1886209374839743,
            1223329526802818,
            105300633354275,
        ]),
        Z: FieldElement51::from_limbs([1, 0, 0, 0, 0]),
        T: FieldElement51::from_limbs([
            803472979097708,
            393902981724766,
            1158077081819914,
            574391322974006,
            336294660666841,
        ]),
    },
    EdwardsPoint {
        X: FieldElement51::from_limbs([0, 0, 0, 0, 0]),
        Y: FieldElement51::from_limbs([
            2251799813685228,
            2251799813685247,
            2251799813685247,
            2251799813685247,
            2251799813685247,
        ]),
        Z: FieldElement51::from_limbs([1, 0, 0, 0, 0]),
        T: FieldElement51::from_limbs([0, 0, 0, 0, 0]),
    },
    EdwardsPoint {
        X: FieldElement51::from_limbs([
            1893055065632419,
            560215195444267,
            1274149604399886,
            821933901047523,
            1691754969406571,
        ]),
        Y: FieldElement51::from_limbs([
            2166873539340326,
            1778179147085316,
            1886209374839743,
            1223329526802818,
            105300633354275,
        ]),
        Z: FieldElement51::from_limbs([1, 0, 0, 0, 0]),
        T: FieldElement51::from_limbs([
            1448326834587521,
            1857896831960481,
            1093722731865333,
            1677408490711241,
            1915505153018406,
        ]),
    },
    EdwardsPoint {
        X: FieldElement51::from_limbs([
            1718705420411056,
            234908883556509,
            2233514472574048,
            2117202627021982,
            765476049583133,
        ]),
        Y: FieldElement51::from_limbs([0, 0, 0, 0, 0]),
        Z: FieldElement51::from_limbs([1, 0, 0, 0, 0]),
        T: FieldElement51::from_limbs([0, 0, 0, 0, 0]),
    },
    EdwardsPoint {
        X: FieldElement51::from_limbs([
            1893055065632419,
            560215195444267,
            1274149604399886,
            821933901047523,
            1691754969406571,
        ]),
        Y: FieldElement51::from_limbs([
            84926274344903,
            473620666599931,
            365590438845504,
            1028470286882429,
            2146499180330972,
        ]),
        Z: FieldElement51::from_limbs([1, 0, 0, 0, 0]),
        T: FieldElement51::from_limbs([
            803472979097708,
            393902981724766,
            1158077081819914,
            574391322974006,
            336294660666841,
        ]),
    },
];

/// Odd multiples of the basepoint `[B, 3B, 5B, 7B, 9B, 11B, 13B, 15B, ..., 127B]`.
#[cfg(feature = "precomputed-tables")]
#[allow(dead_code)]
pub(crate) const AFFINE_ODD_MULTIPLES_OF_BASEPOINT: NafLookupTable8<AffineNielsPoint> =
    NafLookupTable8([
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3540182452943730,
                2497478415033846,
                2521227595762870,
                1462984067271729,
                2389212253076811,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                62697248952638,
                204681361388450,
                631292143396476,
                338455783676468,
                1213667448819585,
            ]),
            xy2d: FieldElement51::from_limbs([
                301289933810280,
                1259582250014073,
                1422107436869536,
                796239922652654,
                1953934009299142,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1601611775252272,
                1720807796594148,
                1132070835939856,
                3512254832574799,
                2147779492816910,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                316559037616741,
                2177824224946892,
                1459442586438991,
                1461528397712656,
                751590696113597,
            ]),
            xy2d: FieldElement51::from_limbs([
                1850748884277385,
                1200145853858453,
                1068094770532492,
                672251375690438,
                1586055907191707,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                769950342298400,
                2384754244604994,
                3095885746880802,
                3225892188161580,
                2977876099231263,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                425251763115706,
                608463272472562,
                442562545713235,
                837766094556764,
                374555092627893,
            ]),
            xy2d: FieldElement51::from_limbs([
                1086255230780037,
                274979815921559,
                1960002765731872,
                929474102396301,
                1190409889297339,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2916800678241215,
                2065379846933858,
                2622030924071124,
                2602788184473875,
                1233371373142984,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                2019367628972465,
                676711900706637,
                110710997811333,
                1108646842542025,
                517791959672113,
            ]),
            xy2d: FieldElement51::from_limbs([
                965130719900578,
                247011430587952,
                526356006571389,
                91986625355052,
                2157223321444601,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1802695059464988,
                1664899123557221,
                2845359304426105,
                2160434469266658,
                3179370264440279,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1725674970513508,
                1933645953859181,
                1542344539275782,
                1767788773573747,
                1297447965928905,
            ]),
            xy2d: FieldElement51::from_limbs([
                1381809363726107,
                1430341051343062,
                2061843536018959,
                1551778050872521,
                2036394857967624,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                4222693909998302,
                2779866139518454,
                1619374932191226,
                2207306624415883,
                1169170329061080,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                2070390218572616,
                1458919061857835,
                624171843017421,
                1055332792707765,
                433987520732508,
            ]),
            xy2d: FieldElement51::from_limbs([
                893653801273833,
                1168026499324677,
                1242553501121234,
                1306366254304474,
                1086752658510815,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2465253816303469,
                3191571337672685,
                1159882208056013,
                2569188183312765,
                621213314200686,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1971678598905747,
                338026507889165,
                762398079972271,
                655096486107477,
                42299032696322,
            ]),
            xy2d: FieldElement51::from_limbs([
                177130678690680,
                1754759263300204,
                1864311296286618,
                1180675631479880,
                1292726903152791,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1913163449625248,
                2712579013977241,
                2193883288642313,
                1008900146920800,
                1721983679009502,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1070401523076875,
                1272492007800961,
                1910153608563310,
                2075579521696771,
                1191169788841221,
            ]),
            xy2d: FieldElement51::from_limbs([
                692896803108118,
                500174642072499,
                2068223309439677,
                1162190621851337,
                1426986007309901,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1819621230288238,
                2735700366193240,
                1755134670739586,
                3080648199451191,
                4172807995775876,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                992069868904071,
                799011518185730,
                1777586403832768,
                1134820506145684,
                1999461475558530,
            ]),
            xy2d: FieldElement51::from_limbs([
                425204543703124,
                2040469794090382,
                1651690622153809,
                1500530168597569,
                1253908377065966,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2105824306960939,
                1387520302709358,
                3633176580451016,
                2211816663841753,
                1629085891776489,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1485201376284999,
                1022406647424656,
                504181009209019,
                962621520820995,
                590876713147230,
            ]),
            xy2d: FieldElement51::from_limbs([
                265873406365287,
                1192742653492898,
                88553098803050,
                525037770869640,
                1266933811251234,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3552316659826612,
                1254279525791875,
                1609927932077699,
                3578654071679972,
                3750681296069893,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                37186803519861,
                1404297334376301,
                578519728836650,
                1740727951192592,
                2095534282477028,
            ]),
            xy2d: FieldElement51::from_limbs([
                833234263154399,
                2023862470013762,
                1854137933982069,
                853924318090959,
                1589812702805850,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3679150557957763,
                1319179453661745,
                497496853611112,
                2665464286942351,
                1208137952365560,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1654513078530905,
                907489875842908,
                126098711296368,
                1726320004173677,
                28269495058173,
            ]),
            xy2d: FieldElement51::from_limbs([
                114436686957443,
                532739313025996,
                115428841215897,
                2191499400074366,
                370280402676434,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1111146849833253,
                2016430049079759,
                1860522747477948,
                3537164738290194,
                4137142824844184,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                429069864577128,
                975327637149449,
                237881983565075,
                1654761232378630,
                2122527599091807,
            ]),
            xy2d: FieldElement51::from_limbs([
                2093793463548278,
                754827233241879,
                1420389751719629,
                1829952782588138,
                2011865756773717,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                676293365438898,
                2850296017886344,
                1205350322490195,
                2763699392265669,
                2133931188538142,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                48340340349120,
                1299261101494832,
                1137329686775218,
                1534848106674340,
                1351662218216799,
            ]),
            xy2d: FieldElement51::from_limbs([
                1904520614137939,
                1590301001714014,
                215781420985270,
                2043534301034629,
                1970888949300424,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2365217962409710,
                2061307169694064,
                1887478590157603,
                2169639621284316,
                2373810867477200,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1020052624656948,
                1260412094216707,
                366721640607121,
                585331442306596,
                345876457758061,
            ]),
            xy2d: FieldElement51::from_limbs([
                975390299880933,
                1066555195234642,
                12651997758352,
                1184252205433068,
                1058378155074223,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1431537716602643,
                2024827957433813,
                3746434518400495,
                1087794891033550,
                2156817571680455,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                929288033346881,
                255179964546973,
                711057989588035,
                208899572612840,
                185348357387383,
            ]),
            xy2d: FieldElement51::from_limbs([
                823689746424808,
                47266130989546,
                209403309368097,
                1100966895202707,
                710792075292719,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2311213117823762,
                3296668540922318,
                2004276520649823,
                1861500579441125,
                3148029033359833,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1563693677475261,
                1843782073741194,
                1950700654453170,
                911540858113949,
                2085151496302359,
            ]),
            xy2d: FieldElement51::from_limbs([
                1427880892005482,
                106216431121745,
                42608394782284,
                1217295886989793,
                1514235272796882,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3544335535746750,
                2367994491347456,
                2567261456502612,
                1854058085060971,
                2263545563461076,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                787426011300053,
                2105981035769060,
                1130476291127206,
                1748659348100075,
                53470983013756,
            ]),
            xy2d: FieldElement51::from_limbs([
                553548273865386,
                5927805718390,
                65184587381926,
                633576679686953,
                576048559439973,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                993787326657446,
                3868807161609258,
                1615796046728943,
                2514644292681953,
                2059021068660907,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                251010270518880,
                1681684095763484,
                1521949356387564,
                431593457045116,
                1855308922422910,
            ]),
            xy2d: FieldElement51::from_limbs([
                618490909691959,
                1257497595618257,
                202952467594088,
                35577762721238,
                1494883566841973,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1673474571932262,
                2409784519770613,
                2636095316260487,
                2761112584601925,
                3333713288149876,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1600640202645197,
                1019569075331823,
                1041916487915822,
                1680448171313267,
                2126903137527901,
            ]),
            xy2d: FieldElement51::from_limbs([
                894964745143659,
                106116880092678,
                1009869382959477,
                317866368542032,
                1986983122763912,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1765281781276487,
                2863247187455184,
                2589075472439062,
                1386435905543054,
                2182338478845320,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1144730936996693,
                2213315231278180,
                1489676672185125,
                665039429138074,
                1131283313040268,
            ]),
            xy2d: FieldElement51::from_limbs([
                2004734176670602,
                1738311085075235,
                418866995976618,
                1050782508034394,
                577747313404652,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2185209688340293,
                1309276076461009,
                2514740038571278,
                3994889904012999,
                3018098826231021,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1405936970888515,
                1754621155316654,
                1211862168554999,
                1813045702919083,
                997853418197172,
            ]),
            xy2d: FieldElement51::from_limbs([
                82037622045021,
                1646398333621944,
                613095452763466,
                1312329542583705,
                81014679202721,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2389287991277873,
                403851022333257,
                1597473361477193,
                2953351602509212,
                2135174663049062,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1826548187201150,
                302299893734126,
                1475477168615781,
                842617616347376,
                1438600873676130,
            ]),
            xy2d: FieldElement51::from_limbs([
                663049852468609,
                1649295727846569,
                1048009692742781,
                628866177992421,
                1914360327429204,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1795645928096646,
                306878154408959,
                2924901319092394,
                2801261341654799,
                1653782432983523,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                2077597317438627,
                212642017882064,
                674844477518888,
                875487498687554,
                2060550250171182,
            ]),
            xy2d: FieldElement51::from_limbs([
                1420448018683809,
                1032663994771382,
                1341927003385267,
                1340360916546159,
                1988547473895228,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1082660122598844,
                2545055705583789,
                3888919679589007,
                1670283344995811,
                3403239134794618,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                90430593339788,
                1838338032241275,
                571293238480915,
                1639938867416883,
                257378872001111,
            ]),
            xy2d: FieldElement51::from_limbs([
                1528535658865034,
                1516636853043960,
                787000569996728,
                1464531394704506,
                1684822625133795,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                811329918113934,
                2783463529007378,
                1769095754634835,
                2970819621866866,
                881037178164325,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1784566501964517,
                433890943689325,
                1186055625589419,
                1496077405487512,
                1731807117886548,
            ]),
            xy2d: FieldElement51::from_limbs([
                424909811816304,
                1355993963741797,
                409606483251841,
                455665350637068,
                1617009023642808,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2478728492077816,
                2780289048655501,
                2328687177473769,
                4107341333582032,
                1316147724308250,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1617420574301156,
                1741273341070467,
                667135503486508,
                2100436564640123,
                1032223920000865,
            ]),
            xy2d: FieldElement51::from_limbs([
                1753947659404033,
                247279202390193,
                1819288880178945,
                737334285670249,
                1037873664856104,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1762568490530034,
                673742465299012,
                2054571050635888,
                2040165159255111,
                3040123733327257,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1627187989987422,
                1686331580821752,
                1309895873498183,
                719718719104086,
                300063199808722,
            ]),
            xy2d: FieldElement51::from_limbs([
                238176707016164,
                1440454788877048,
                203336037573144,
                1437789888677072,
                101522256664211,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1895216760098480,
                1934324337975022,
                3677350688973167,
                2536415965456176,
                714678003308640,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                508185358728815,
                1691320535341855,
                2168887448239256,
                1035124393070661,
                1936603999698584,
            ]),
            xy2d: FieldElement51::from_limbs([
                390562831571647,
                1390223890708972,
                1383183990676371,
                435998174196410,
                1882086414390730,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3747620842612921,
                2081794785291195,
                3284594056262745,
                2090090346797895,
                2581692978935809,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                244144781251265,
                1290834426417077,
                1888701171101942,
                1233922456644870,
                241117402207491,
            ]),
            xy2d: FieldElement51::from_limbs([
                1266169390045455,
                1148042013187970,
                878921907853942,
                1815738019658093,
                908920199341621,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2521768507305118,
                953557056811112,
                2015863732865770,
                1358382511861315,
                2835421647899992,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                2239837206240498,
                330928973149665,
                422268062913642,
                1481280019493032,
                619879520439841,
            ]),
            xy2d: FieldElement51::from_limbs([
                1360166735366017,
                1770556573948510,
                1395061284191031,
                1814003148068126,
                522781147076884,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2611794802645686,
                707234844948070,
                1314059396506491,
                2919250341703934,
                2161831667832785,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                934831784182383,
                433734253968318,
                1660867106725771,
                1968393082772831,
                873946300968490,
            ]),
            xy2d: FieldElement51::from_limbs([
                26306827827554,
                430884999378685,
                1504310424376419,
                1761358720837522,
                542195685418530,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1762131062631725,
                3123952634417535,
                3619918390837537,
                2909990877347294,
                1411594230004385,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                538272372224622,
                1425714779586199,
                588313661410172,
                1497062084392578,
                1602174047128512,
            ]),
            xy2d: FieldElement51::from_limbs([
                907490361939255,
                1963620338391363,
                626927432296975,
                1250748516081414,
                959901171882527,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1335066153744413,
                2887804660779657,
                2653073855954038,
                2765226981667422,
                938831784476763,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                296699434737224,
                2047543711075683,
                2076451038937139,
                227783599906901,
                1602062110967627,
            ]),
            xy2d: FieldElement51::from_limbs([
                1574834773194203,
                1384279952062839,
                393652417255803,
                2166968242848859,
                1552890441390820,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1619646774410947,
                1576090644023562,
                3035228391320965,
                1735328519940543,
                2355324535937066,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1024074573633446,
                957088456885874,
                1690425531356997,
                2102187380180052,
                1082544623222033,
            ]),
            xy2d: FieldElement51::from_limbs([
                1871906170635853,
                1719383891167200,
                1584032250247862,
                823764804192117,
                2244048510084261,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                642147846489775,
                3334304977145699,
                305205716788147,
                2589176626729533,
                2224680511484174,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1734162377166545,
                260713621840346,
                157174591942595,
                952544272517991,
                222818702471733,
            ]),
            xy2d: FieldElement51::from_limbs([
                1213115494182947,
                286778704335711,
                2130189536016490,
                308349182281342,
                1217623948685491,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3360052266973635,
                1843486583624091,
                1561693837124349,
                1084041964025479,
                1866270922024009,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                460705465481210,
                1968151453817859,
                497005926994844,
                625618055866751,
                2176893440866887,
            ]),
            xy2d: FieldElement51::from_limbs([
                1655800250476757,
                2036588542300609,
                666447448675243,
                1615721995750683,
                1508669225186765,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2245948203759141,
                1058306669699396,
                1452898014240582,
                3961024141962768,
                1633235287338608,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                986647273684279,
                1507266907811370,
                1260572633649005,
                2071672342077446,
                695976026010857,
            ]),
            xy2d: FieldElement51::from_limbs([
                1312356620823495,
                1635278548098567,
                901946076841033,
                585120475533168,
                1240667113237384,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2313723935779695,
                1506054666773895,
                996040223525031,
                636592914999692,
                1497801917020297,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                292042016419794,
                1158932298133044,
                2062611870323738,
                1946058478962569,
                1749165808126286,
            ]),
            xy2d: FieldElement51::from_limbs([
                654683942212830,
                1526897351349087,
                2006818439922838,
                2194919327350361,
                1451960776874416,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3015041017808905,
                2951823141773809,
                2584865668253675,
                2508192032998563,
                2582137700042019,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1628123495344283,
                2072923641214546,
                1647225812023982,
                855655925244679,
                1758126430071140,
            ]),
            xy2d: FieldElement51::from_limbs([
                1615895096489599,
                275295258643784,
                937665541219916,
                1313496726746346,
                1186468946422626,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1603070202850694,
                2072127623773242,
                1692648737212158,
                2493373404187852,
                1248948672117105,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                11167836031898,
                596565174397990,
                2196351068723859,
                314744641791907,
                1102014997250781,
            ]),
            xy2d: FieldElement51::from_limbs([
                1409047922401191,
                69960384467966,
                688103515547600,
                1309746102488044,
                150292892873778,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1986083055103168,
                691715819340300,
                1361811659746933,
                3459052030333434,
                1063594696046061,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1201987338414749,
                2198784582460616,
                1203335513981498,
                489243077045066,
                2205278143582433,
            ]),
            xy2d: FieldElement51::from_limbs([
                2034744376624534,
                2077387101466387,
                148448542974969,
                1502697574577258,
                473186584705655,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                472016956315960,
                720786972252993,
                2840633661190043,
                3150798753357827,
                2816563335499153,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                253464247569755,
                168314237403057,
                511780806170295,
                1058862316549135,
                1646858476817137,
            ]),
            xy2d: FieldElement51::from_limbs([
                595092995922219,
                1491311840717691,
                291581784452778,
                1569186646367854,
                1031385061400544,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3483137021572755,
                1526955102024322,
                2778006642704458,
                457549634924205,
                1097420237736736,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1246991699537710,
                81367319519439,
                530844036072196,
                163656863755855,
                1950742455979290,
            ]),
            xy2d: FieldElement51::from_limbs([
                191532664076407,
                539378506082089,
                1021612562876554,
                1026603384732632,
                1773368780410653,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                4144620731387879,
                590179521333342,
                4034023318016108,
                2255745030335426,
                2699746851701250,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                2206599697359952,
                553895797384417,
                181689161933786,
                1153123447919104,
                778568064152659,
            ]),
            xy2d: FieldElement51::from_limbs([
                1706307000059211,
                1885601289314487,
                889758608505788,
                550131729999853,
                1006862664714268,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3210197754285058,
                2048500453422630,
                3403309827888207,
                927154428508963,
                4199813798872019,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                992058915374933,
                476120535358775,
                1973648780784340,
                2025282643598818,
                2182318983793230,
            ]),
            xy2d: FieldElement51::from_limbs([
                1343440812005821,
                1316045839091795,
                1884951299078063,
                1765919609219175,
                2197567554627988,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3129247779382818,
                4415026969054274,
                1900265885969643,
                1528796215447059,
                2172730393748688,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1773355092297603,
                64654329538271,
                1332124041660957,
                748492100858001,
                895500006200535,
            ]),
            xy2d: FieldElement51::from_limbs([
                2000840647851980,
                546565968824914,
                420633283457524,
                195470736374507,
                1958689297569520,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                743138980705446,
                3411117504637167,
                2591389959690621,
                2380042066577202,
                3022267940115114,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                165947002229363,
                115186103724967,
                1068573292121517,
                1842565776920938,
                1969395681111987,
            ]),
            xy2d: FieldElement51::from_limbs([
                553322266190633,
                234265665613185,
                484544650202821,
                1238773526575826,
                2017991917953668,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2581954631514051,
                1245093644265357,
                3537016673825374,
                1834216551713857,
                923978372152807,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1855378315339552,
                890045579230758,
                1764718173975590,
                197904186055854,
                1718129022310327,
            ]),
            xy2d: FieldElement51::from_limbs([
                1278162928734862,
                1894118254109862,
                987503995465517,
                177406744098996,
                781538103127693,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1996603431230215,
                1191888797552937,
                1207440075928499,
                2765853449051137,
                2525314961343288,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                808903879370889,
                990820108751280,
                1084429472258867,
                1078562781312589,
                254514692695625,
            ]),
            xy2d: FieldElement51::from_limbs([
                615855140068469,
                586046731175395,
                693470779212674,
                1964537100203868,
                1350330550265229,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3344544372023708,
                720386671449874,
                2480841360702110,
                2036034126860286,
                2015744690201389,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1337446193390478,
                1984110761311871,
                746489405020285,
                407347127604128,
                1740475330360596,
            ]),
            xy2d: FieldElement51::from_limbs([
                140840424783613,
                1063284623568331,
                1136446106453878,
                372042229029799,
                442607248430694,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2330781679120937,
                376801425148230,
                2032603686676107,
                1488926293635130,
                1317278311532959,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1290116731380016,
                2166899563471713,
                831997001838078,
                870954980505220,
                2108537278055823,
            ]),
            xy2d: FieldElement51::from_limbs([
                1912719171026343,
                846194720551034,
                2043988124740726,
                993234269653961,
                421229796383281,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2651184584992902,
                2775702557638963,
                2539786009779572,
                2575974880015305,
                2122619079836732,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1154054290132562,
                931753998725577,
                1647742001778052,
                865765466488226,
                1083816107290025,
            ]),
            xy2d: FieldElement51::from_limbs([
                986341121095108,
                1522330369638573,
                1990880546211047,
                501525962272123,
                198539304862139,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1496414019192687,
                3991034436173951,
                3380311659062196,
                2854747485359158,
                3346958036643152,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                805612068303425,
                1891790027761335,
                1587008567571549,
                722120737390201,
                378156757163816,
            ]),
            xy2d: FieldElement51::from_limbs([
                1588994517921951,
                977362751042302,
                1329302387067714,
                2069348224564088,
                1586007159625211,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2490539421551682,
                1985699850375015,
                2331762317128172,
                4145097393776678,
                2521049460190674,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                615817553313996,
                2245962768078178,
                482564324326173,
                2101336843140780,
                1240914880829407,
            ]),
            xy2d: FieldElement51::from_limbs([
                1438242482238189,
                874267817785463,
                1620810389770625,
                866155221338671,
                1040426546798301,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                2403083624110300,
                2548561409802975,
                2492699136535911,
                2358289519456539,
                3203964320363148,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1913986535403097,
                1977163223054199,
                1972905914623196,
                1650122133472502,
                1905849310819035,
            ]),
            xy2d: FieldElement51::from_limbs([
                858174816360838,
                614595356564037,
                1099584959044836,
                636998087084906,
                1070393269058348,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3666695924830668,
                3585640662737501,
                2372994528684236,
                2628565977288995,
                3482812783469694,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1994161359147952,
                2198039369802658,
                62790022842537,
                1522306785848169,
                951223194802833,
            ]),
            xy2d: FieldElement51::from_limbs([
                852296621440717,
                431889737774209,
                370755457746189,
                437604073958073,
                627857326892757,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1794955764684156,
                2586904290013612,
                1322647643615887,
                856117964085888,
                2652432778663153,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                933592377399646,
                78031722952813,
                926049890685253,
                1471649501316246,
                33789909190376,
            ]),
            xy2d: FieldElement51::from_limbs([
                1479319468832059,
                203906207621608,
                659828362330083,
                44358398435755,
                1273573524210803,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1592342143350813,
                3227219208247713,
                2345240352078765,
                2577750109932929,
                2933512841197243,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                2184946892642995,
                1517382324576002,
                1557940277419806,
                2170635134813213,
                747314658627002,
            ]),
            xy2d: FieldElement51::from_limbs([
                1823193620577742,
                1135817878516419,
                1731253819308581,
                1031652967267804,
                2123506616999453,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1346190246005805,
                2052692552023851,
                1718128041785940,
                2491557332978474,
                3474370880388305,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                424776012994573,
                281050757243423,
                626466040846420,
                990194703866532,
                38571969885982,
            ]),
            xy2d: FieldElement51::from_limbs([
                192408346595466,
                1054889725292349,
                584097975693004,
                1447909807397749,
                2134645004369136,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3169895788615063,
                3503097743181446,
                601598510029975,
                1422812237223371,
                2121009661378329,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1603348391996783,
                2066143816131699,
                1789627290363958,
                2145705961178118,
                1985578641438222,
            ]),
            xy2d: FieldElement51::from_limbs([
                352633958653380,
                856927627345554,
                793925083122702,
                93551575767286,
                1222010153634215,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1756866499986349,
                911731956999969,
                2707505543214075,
                4006920335263786,
                822501008147910,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1094036422864347,
                1897208881572508,
                1503607738246960,
                1901060196071406,
                294068411105729,
            ]),
            xy2d: FieldElement51::from_limbs([
                587776484399576,
                1116861711228807,
                343398777436088,
                936544065763093,
                1643746750211060,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                3477749685790410,
                267997399528836,
                2953780922004404,
                3252368924080907,
                3787792887348381,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                2042368155872443,
                41662387210459,
                1676313264498480,
                1333968523426810,
                1765708383352310,
            ]),
            xy2d: FieldElement51::from_limbs([
                1453394896690938,
                1585795827439909,
                1469309456804303,
                1294645324464404,
                2042954198665899,
            ]),
        },
        AffineNielsPoint {
            y_plus_x: FieldElement51::from_limbs([
                1810069207599881,
                1358344669503239,
                1989371257548167,
                2316270051121225,
                3019675451276507,
            ]),
            y_minus_x: FieldElement51::from_limbs([
                1866114438287676,
                1663420339568364,
                1437691317033088,
                538298302628038,
                1212711449614363,
            ]),
            xy2d: FieldElement51::from_limbs([
                1769235035677897,
                1562012115317882,
                31277513664750,
                536198657928416,
                1976134212537183,
            ]),
        },
    ]);
