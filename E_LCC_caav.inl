namespace MRLCC_CAAV {

  FTensorDecl TensorDecls[37] = {
    /*  0*/{"R"      , "cAae"     , "", USAGE_Residual    },
    /*  1*/{"t"      , "cAae"     , "", USAGE_Amplitude   },
    /*  2*/{"T"      , "cAae"     , "", USAGE_Amplitude   },
    /*  3*/{"b"      , "cAae"     , "", USAGE_Amplitude   },
    /*  4*/{"B"      , "cAae"     , "", USAGE_Amplitude   },
    /*  5*/{"p"      , "cAae"     , "", USAGE_Amplitude   },
    /*  6*/{"Ap"     , "cAae"     , "", USAGE_Amplitude   },
    /*  7*/{"P"      , "cAae"     , "", USAGE_Amplitude   },
    /*  8*/{"AP"     , "cAae"     , "", USAGE_Amplitude   },
    /*  9*/{"p1"     , "caae"     , "", USAGE_Amplitude   },
    /* 10*/{"p2"     , "caae"     , "", USAGE_Amplitude   },
    /* 11*/{"Ap1"    , "caae"     , "", USAGE_Amplitude   },
    /* 12*/{"Ap2"    , "caae"     , "", USAGE_Amplitude   },
    /* 13*/{"b1"     , "caae"     , "", USAGE_Amplitude   },
    /* 14*/{"b2"     , "caae"     , "", USAGE_Amplitude   },
    /* 15*/{"B1"     , "caae"     , "", USAGE_Amplitude   },
    /* 16*/{"B2"     , "caae"     , "", USAGE_Amplitude   },
    /* 17*/{"f"      , "ec"       , "", USAGE_Hamiltonian },
    /* 18*/{"k"      , "aa"       , "", USAGE_Hamiltonian },
    /* 19*/{"k"      , "cc"       , "", USAGE_Hamiltonian },
    /* 20*/{"k"      , "ee"       , "", USAGE_Hamiltonian },
    /* 21*/{"W"      , "aaaa"     , "", USAGE_Hamiltonian },
    /* 22*/{"W"      , "aeae"     , "", USAGE_Hamiltonian },
    /* 23*/{"W"      , "aeea"     , "", USAGE_Hamiltonian },
    /* 24*/{"W"      , "cccc"     , "", USAGE_Hamiltonian },
    /* 25*/{"W"      , "caca"     , "", USAGE_Hamiltonian },
    /* 26*/{"W"      , "caac"     , "", USAGE_Hamiltonian },
    /* 27*/{"W"      , "cece"     , "", USAGE_Hamiltonian },
    /* 28*/{"W"      , "ceec"     , "", USAGE_Hamiltonian },
    /* 29*/{"W"      , "eaca"     , "", USAGE_Hamiltonian },
    /* 30*/{"W"      , "aeca"     , "", USAGE_Hamiltonian },
    /* 31*/{"delta"  , "aa"       , "", USAGE_Density     },
    /* 32*/{"delta"  , "cc"       , "", USAGE_Density     },
    /* 33*/{"E1"     , "aa"       , "", USAGE_Density     },
    /* 34*/{"E2"     , "aaaa"     , "", USAGE_Density     },
    /* 35*/{"S1"     , "AaAa"     , "", USAGE_Density     },
    /* 36*/{"E3"     , "aaaaaa"   , "", USAGE_Density     },
  };

	FEqInfo EqsRes[216] = {
    {"abcd,be,aefd,cf",          2.0,   4, {11,18, 9,33}},     //   2.0 Ap1[abcd] k[be] p1[aefd] E1[cf] 
    {"abcd,ef,abed,cf",         -2.0,   4, {11,18, 9,33}},     //  -2.0 Ap1[abcd] k[ef] p1[abed] E1[cf] 
    {"abcd,be,aefd,cf",         -1.0,   4, {11,18,10,33}},     //  -1.0 Ap1[abcd] k[be] p2[aefd] E1[cf] 
    {"abcd,ef,abed,cf",          1.0,   4, {11,18,10,33}},     //   1.0 Ap1[abcd] k[ef] p2[abed] E1[cf] 
    {"abcd,ae,ebfd,cf",         -2.0,   4, {11,19, 9,33}},     //  -2.0 Ap1[abcd] k[ae] p1[ebfd] E1[cf] 
    {"abcd,ae,ebfd,cf",          1.0,   4, {11,19,10,33}},     //   1.0 Ap1[abcd] k[ae] p2[ebfd] E1[cf] 
    {"abcd,de,abfe,cf",          2.0,   4, {11,20, 9,33}},     //   2.0 Ap1[abcd] k[de] p1[abfe] E1[cf] 
    {"abcd,de,abfe,cf",         -1.0,   4, {11,20,10,33}},     //  -1.0 Ap1[abcd] k[de] p2[abfe] E1[cf] 
    {"abcd,bdef,aegf,cg",        2.0,   4, {11,22, 9,33}},     //   2.0 Ap1[abcd] W[bdef] p1[aegf] E1[cg] 
    {"abcd,bdef,aegf,cg",       -1.0,   4, {11,22,10,33}},     //  -1.0 Ap1[abcd] W[bdef] p2[aegf] E1[cg] 
    {"abcd,bdef,afge,cg",       -1.0,   4, {11,23, 9,33}},     //  -1.0 Ap1[abcd] W[bdef] p1[afge] E1[cg] 
    {"abcd,bdef,afge,cg",        2.0,   4, {11,23,10,33}},     //   2.0 Ap1[abcd] W[bdef] p2[afge] E1[cg] 
    {"abcd,aehf,fbgd,cg,he",     2.0,   5, {11,24, 9,33,32}},  //   2.0 Ap1[abcd] W[aehf] p1[fbgd] E1[cg] delta[he] 
    {"abcd,aefh,fbgd,cg,he",    -4.0,   5, {11,24, 9,33,32}},  //  -4.0 Ap1[abcd] W[aefh] p1[fbgd] E1[cg] delta[he] 
    {"abcd,aehf,fbgd,cg,he",    -1.0,   5, {11,24,10,33,32}},  //  -1.0 Ap1[abcd] W[aehf] p2[fbgd] E1[cg] delta[he] 
    {"abcd,aefh,fbgd,cg,he",     2.0,   5, {11,24,10,33,32}},  //   2.0 Ap1[abcd] W[aefh] p2[fbgd] E1[cg] delta[he] 
    {"abcd,abef,efgd,cg",       -2.0,   4, {11,25, 9,33}},     //  -2.0 Ap1[abcd] W[abef] p1[efgd] E1[cg] 
    {"abcd,ebhf,afgd,cg,he",     4.0,   5, {11,25, 9,33,32}},  //   4.0 Ap1[abcd] W[ebhf] p1[afgd] E1[cg] delta[he] 
    {"abcd,efhg,abfd,cg,he",    -4.0,   5, {11,25, 9,33,32}},  //  -4.0 Ap1[abcd] W[efhg] p1[abfd] E1[cg] delta[he] 
    {"abcd,abef,efgd,cg",        1.0,   4, {11,25,10,33}},     //   1.0 Ap1[abcd] W[abef] p2[efgd] E1[cg] 
    {"abcd,ebhf,afgd,cg,he",    -2.0,   5, {11,25,10,33,32}},  //  -2.0 Ap1[abcd] W[ebhf] p2[afgd] E1[cg] delta[he] 
    {"abcd,efhg,abfd,cg,he",     2.0,   5, {11,25,10,33,32}},  //   2.0 Ap1[abcd] W[efhg] p2[abfd] E1[cg] delta[he] 
    {"abcd,ebfa,efgd,cg",        1.0,   4, {11,26, 9,33}},     //   1.0 Ap1[abcd] W[ebfa] p1[efgd] E1[cg] 
    {"abcd,ebfh,afgd,cg,he",    -2.0,   5, {11,26, 9,33,32}},  //  -2.0 Ap1[abcd] W[ebfh] p1[afgd] E1[cg] delta[he] 
    {"abcd,efgh,abfd,cg,he",     2.0,   5, {11,26, 9,33,32}},  //   2.0 Ap1[abcd] W[efgh] p1[abfd] E1[cg] delta[he] 
    {"abcd,ebfa,efgd,cg",       -2.0,   4, {11,26,10,33}},     //  -2.0 Ap1[abcd] W[ebfa] p2[efgd] E1[cg] 
    {"abcd,ebfh,afgd,cg,he",     1.0,   5, {11,26,10,33,32}},  //   1.0 Ap1[abcd] W[ebfh] p2[afgd] E1[cg] delta[he] 
    {"abcd,efgh,abfd,cg,he",    -1.0,   5, {11,26,10,33,32}},  //  -1.0 Ap1[abcd] W[efgh] p2[abfd] E1[cg] delta[he] 
    {"abcd,adef,ebgf,cg",       -2.0,   4, {11,27, 9,33}},     //  -2.0 Ap1[abcd] W[adef] p1[ebgf] E1[cg] 
    {"abcd,edhf,abgf,cg,he",     4.0,   5, {11,27, 9,33,32}},  //   4.0 Ap1[abcd] W[edhf] p1[abgf] E1[cg] delta[he] 
    {"abcd,adef,ebgf,cg",        1.0,   4, {11,27,10,33}},     //   1.0 Ap1[abcd] W[adef] p2[ebgf] E1[cg] 
    {"abcd,edhf,abgf,cg,he",    -2.0,   5, {11,27,10,33,32}},  //  -2.0 Ap1[abcd] W[edhf] p2[abgf] E1[cg] delta[he] 
    {"abcd,edfa,ebgf,cg",        4.0,   4, {11,28, 9,33}},     //   4.0 Ap1[abcd] W[edfa] p1[ebgf] E1[cg] 
    {"abcd,edfh,abgf,cg,he",    -2.0,   5, {11,28, 9,33,32}},  //  -2.0 Ap1[abcd] W[edfh] p1[abgf] E1[cg] delta[he] 
    {"abcd,edfa,ebgf,cg",       -2.0,   4, {11,28,10,33}},     //  -2.0 Ap1[abcd] W[edfa] p2[ebgf] E1[cg] 
    {"abcd,edfh,abgf,cg,he",     1.0,   5, {11,28,10,33,32}},  //   1.0 Ap1[abcd] W[edfh] p2[abgf] E1[cg] delta[he] 
    {"abcd,be,aefd,cf",         -1.0,   4, {12,18, 9,33}},     //  -1.0 Ap2[abcd] k[be] p1[aefd] E1[cf] 
    {"abcd,ef,abed,cf",          1.0,   4, {12,18, 9,33}},     //   1.0 Ap2[abcd] k[ef] p1[abed] E1[cf] 
    {"abcd,be,aefd,cf",          2.0,   4, {12,18,10,33}},     //   2.0 Ap2[abcd] k[be] p2[aefd] E1[cf] 
    {"abcd,ef,abed,cf",         -2.0,   4, {12,18,10,33}},     //  -2.0 Ap2[abcd] k[ef] p2[abed] E1[cf] 
    {"abcd,ae,ebfd,cf",          1.0,   4, {12,19, 9,33}},     //   1.0 Ap2[abcd] k[ae] p1[ebfd] E1[cf] 
    {"abcd,ae,ebfd,cf",         -2.0,   4, {12,19,10,33}},     //  -2.0 Ap2[abcd] k[ae] p2[ebfd] E1[cf] 
    {"abcd,de,abfe,cf",         -1.0,   4, {12,20, 9,33}},     //  -1.0 Ap2[abcd] k[de] p1[abfe] E1[cf] 
    {"abcd,de,abfe,cf",          2.0,   4, {12,20,10,33}},     //   2.0 Ap2[abcd] k[de] p2[abfe] E1[cf] 
    {"abcd,bdef,aegf,cg",       -1.0,   4, {12,22, 9,33}},     //  -1.0 Ap2[abcd] W[bdef] p1[aegf] E1[cg] 
    {"abcd,bdef,aegf,cg",        2.0,   4, {12,22,10,33}},     //   2.0 Ap2[abcd] W[bdef] p2[aegf] E1[cg] 
    {"abcd,bdef,afge,cg",        2.0,   4, {12,23, 9,33}},     //   2.0 Ap2[abcd] W[bdef] p1[afge] E1[cg] 
    {"abcd,bdef,afge,cg",       -1.0,   4, {12,23,10,33}},     //  -1.0 Ap2[abcd] W[bdef] p2[afge] E1[cg] 
    {"abcd,aehf,fbgd,cg,he",    -1.0,   5, {12,24, 9,33,32}},  //  -1.0 Ap2[abcd] W[aehf] p1[fbgd] E1[cg] delta[he] 
    {"abcd,aefh,fbgd,cg,he",     2.0,   5, {12,24, 9,33,32}},  //   2.0 Ap2[abcd] W[aefh] p1[fbgd] E1[cg] delta[he] 
    {"abcd,aehf,fbgd,cg,he",     2.0,   5, {12,24,10,33,32}},  //   2.0 Ap2[abcd] W[aehf] p2[fbgd] E1[cg] delta[he] 
    {"abcd,aefh,fbgd,cg,he",    -4.0,   5, {12,24,10,33,32}},  //  -4.0 Ap2[abcd] W[aefh] p2[fbgd] E1[cg] delta[he] 
    {"abcd,abef,efgd,cg",        1.0,   4, {12,25, 9,33}},     //   1.0 Ap2[abcd] W[abef] p1[efgd] E1[cg] 
    {"abcd,ebhf,afgd,cg,he",    -2.0,   5, {12,25, 9,33,32}},  //  -2.0 Ap2[abcd] W[ebhf] p1[afgd] E1[cg] delta[he] 
    {"abcd,efhg,abfd,cg,he",     2.0,   5, {12,25, 9,33,32}},  //   2.0 Ap2[abcd] W[efhg] p1[abfd] E1[cg] delta[he] 
    {"abcd,abef,efgd,cg",       -2.0,   4, {12,25,10,33}},     //  -2.0 Ap2[abcd] W[abef] p2[efgd] E1[cg] 
    {"abcd,ebhf,afgd,cg,he",     4.0,   5, {12,25,10,33,32}},  //   4.0 Ap2[abcd] W[ebhf] p2[afgd] E1[cg] delta[he] 
    {"abcd,efhg,abfd,cg,he",    -4.0,   5, {12,25,10,33,32}},  //  -4.0 Ap2[abcd] W[efhg] p2[abfd] E1[cg] delta[he] 
    {"abcd,ebfa,efgd,cg",       -2.0,   4, {12,26, 9,33}},     //  -2.0 Ap2[abcd] W[ebfa] p1[efgd] E1[cg] 
    {"abcd,ebfh,afgd,cg,he",     1.0,   5, {12,26, 9,33,32}},  //   1.0 Ap2[abcd] W[ebfh] p1[afgd] E1[cg] delta[he] 
    {"abcd,efgh,abfd,cg,he",    -1.0,   5, {12,26, 9,33,32}},  //  -1.0 Ap2[abcd] W[efgh] p1[abfd] E1[cg] delta[he] 
    {"abcd,ebfa,efgd,cg",        4.0,   4, {12,26,10,33}},     //   4.0 Ap2[abcd] W[ebfa] p2[efgd] E1[cg] 
    {"abcd,ebfh,afgd,cg,he",    -2.0,   5, {12,26,10,33,32}},  //  -2.0 Ap2[abcd] W[ebfh] p2[afgd] E1[cg] delta[he] 
    {"abcd,efgh,abfd,cg,he",     2.0,   5, {12,26,10,33,32}},  //   2.0 Ap2[abcd] W[efgh] p2[abfd] E1[cg] delta[he] 
    {"abcd,adef,ebgf,cg",        1.0,   4, {12,27, 9,33}},     //   1.0 Ap2[abcd] W[adef] p1[ebgf] E1[cg] 
    {"abcd,edhf,abgf,cg,he",    -2.0,   5, {12,27, 9,33,32}},  //  -2.0 Ap2[abcd] W[edhf] p1[abgf] E1[cg] delta[he] 
    {"abcd,adef,ebgf,cg",       -2.0,   4, {12,27,10,33}},     //  -2.0 Ap2[abcd] W[adef] p2[ebgf] E1[cg] 
    {"abcd,edhf,abgf,cg,he",     4.0,   5, {12,27,10,33,32}},  //   4.0 Ap2[abcd] W[edhf] p2[abgf] E1[cg] delta[he] 
    {"abcd,edfa,ebgf,cg",       -2.0,   4, {12,28, 9,33}},     //  -2.0 Ap2[abcd] W[edfa] p1[ebgf] E1[cg] 
    {"abcd,edfh,abgf,cg,he",     1.0,   5, {12,28, 9,33,32}},  //   1.0 Ap2[abcd] W[edfh] p1[abgf] E1[cg] delta[he] 
    {"abcd,edfa,ebgf,cg",        1.0,   4, {12,28,10,33}},     //   1.0 Ap2[abcd] W[edfa] p2[ebgf] E1[cg] 
    {"abcd,edfh,abgf,cg,he",    -2.0,   5, {12,28,10,33,32}},  //  -2.0 Ap2[abcd] W[edfh] p2[abgf] E1[cg] delta[he] 
    {"abcd,ef,aegd,cfbg",        2.0,   4, {11,18, 9,34}},     //   2.0 Ap1[abcd] k[ef] p1[aegd] E2[cfbg] 
    {"abcd,ef,aged,cgbf",       -2.0,   4, {11,18, 9,34}},     //  -2.0 Ap1[abcd] k[ef] p1[aged] E2[cgbf] 
    {"abcd,ef,aegd,cfbg",       -1.0,   4, {11,18,10,34}},     //  -1.0 Ap1[abcd] k[ef] p2[aegd] E2[cfbg] 
    {"abcd,ef,aged,cgbf",        1.0,   4, {11,18,10,34}},     //   1.0 Ap1[abcd] k[ef] p2[aged] E2[cgbf] 
    {"abcd,ae,efgd,cfbg",       -2.0,   4, {11,19, 9,34}},     //  -2.0 Ap1[abcd] k[ae] p1[efgd] E2[cfbg] 
    {"abcd,ae,efgd,cfbg",        1.0,   4, {11,19,10,34}},     //   1.0 Ap1[abcd] k[ae] p2[efgd] E2[cfbg] 
    {"abcd,de,afge,cfbg",        2.0,   4, {11,20, 9,34}},     //   2.0 Ap1[abcd] k[de] p1[afge] E2[cfbg] 
    {"abcd,de,afge,cfbg",       -1.0,   4, {11,20,10,34}},     //  -1.0 Ap1[abcd] k[de] p2[afge] E2[cfbg] 
    {"abcd,befg,aehd,cgfh",      2.0,   4, {11,21, 9,34}},     //   2.0 Ap1[abcd] W[befg] p1[aehd] E2[cgfh] 
    {"abcd,befg,afhd,cehg",      2.0,   4, {11,21, 9,34}},     //   2.0 Ap1[abcd] W[befg] p1[afhd] E2[cehg] 
    {"abcd,befg,ahed,chfg",     -2.0,   4, {11,21, 9,34}},     //  -2.0 Ap1[abcd] W[befg] p1[ahed] E2[chfg] 
    {"abcd,efgh,abed,cfgh",     -2.0,   4, {11,21, 9,34}},     //  -2.0 Ap1[abcd] W[efgh] p1[abed] E2[cfgh] 
    {"abcd,befg,aehd,cgfh",     -1.0,   4, {11,21,10,34}},     //  -1.0 Ap1[abcd] W[befg] p2[aehd] E2[cgfh] 
    {"abcd,befg,afhd,cehg",     -1.0,   4, {11,21,10,34}},     //  -1.0 Ap1[abcd] W[befg] p2[afhd] E2[cehg] 
    {"abcd,befg,ahed,chfg",      1.0,   4, {11,21,10,34}},     //   1.0 Ap1[abcd] W[befg] p2[ahed] E2[chfg] 
    {"abcd,efgh,abed,cfgh",      1.0,   4, {11,21,10,34}},     //   1.0 Ap1[abcd] W[efgh] p2[abed] E2[cfgh] 
    {"abcd,bdef,aghf,cgeh",      2.0,   4, {11,22, 9,34}},     //   2.0 Ap1[abcd] W[bdef] p1[aghf] E2[cgeh] 
    {"abcd,edfg,abhg,cehf",      2.0,   4, {11,22, 9,34}},     //   2.0 Ap1[abcd] W[edfg] p1[abhg] E2[cehf] 
    {"abcd,edfg,aehg,cfbh",      2.0,   4, {11,22, 9,34}},     //   2.0 Ap1[abcd] W[edfg] p1[aehg] E2[cfbh] 
    {"abcd,bdef,aghf,cgeh",     -1.0,   4, {11,22,10,34}},     //  -1.0 Ap1[abcd] W[bdef] p2[aghf] E2[cgeh] 
    {"abcd,edfg,abhg,cehf",     -1.0,   4, {11,22,10,34}},     //  -1.0 Ap1[abcd] W[edfg] p2[abhg] E2[cehf] 
    {"abcd,edfg,aehg,cfbh",     -1.0,   4, {11,22,10,34}},     //  -1.0 Ap1[abcd] W[edfg] p2[aehg] E2[cfbh] 
    {"abcd,bdef,aghe,cgfh",     -1.0,   4, {11,23, 9,34}},     //  -1.0 Ap1[abcd] W[bdef] p1[aghe] E2[cgfh] 
    {"abcd,edfg,abhf,cehg",     -1.0,   4, {11,23, 9,34}},     //  -1.0 Ap1[abcd] W[edfg] p1[abhf] E2[cehg] 
    {"abcd,edfg,aghf,cebh",     -1.0,   4, {11,23, 9,34}},     //  -1.0 Ap1[abcd] W[edfg] p1[aghf] E2[cebh] 
    {"abcd,bdef,aghe,cghf",     -1.0,   4, {11,23,10,34}},     //  -1.0 Ap1[abcd] W[bdef] p2[aghe] E2[cghf] 
    {"abcd,edfg,abhf,cegh",     -1.0,   4, {11,23,10,34}},     //  -1.0 Ap1[abcd] W[edfg] p2[abhf] E2[cegh] 
    {"abcd,edfg,aghf,cebh",      2.0,   4, {11,23,10,34}},     //   2.0 Ap1[abcd] W[edfg] p2[aghf] E2[cebh] 
    {"abcd,aeif,fghd,cgbh,ie",   2.0,   5, {11,24, 9,34,32}},  //   2.0 Ap1[abcd] W[aeif] p1[fghd] E2[cgbh] delta[ie] 
    {"abcd,aefi,fghd,cgbh,ie",  -4.0,   5, {11,24, 9,34,32}},  //  -4.0 Ap1[abcd] W[aefi] p1[fghd] E2[cgbh] delta[ie] 
    {"abcd,aeif,fghd,cgbh,ie",  -1.0,   5, {11,24,10,34,32}},  //  -1.0 Ap1[abcd] W[aeif] p2[fghd] E2[cgbh] delta[ie] 
    {"abcd,aefi,fghd,cgbh,ie",   2.0,   5, {11,24,10,34,32}},  //   2.0 Ap1[abcd] W[aefi] p2[fghd] E2[cgbh] delta[ie] 
    {"abcd,abef,eghd,cgfh",     -2.0,   4, {11,25, 9,34}},     //  -2.0 Ap1[abcd] W[abef] p1[eghd] E2[cgfh] 
    {"abcd,aefg,fbhd,cehg",     -2.0,   4, {11,25, 9,34}},     //  -2.0 Ap1[abcd] W[aefg] p1[fbhd] E2[cehg] 
    {"abcd,aefg,fehd,cgbh",     -2.0,   4, {11,25, 9,34}},     //  -2.0 Ap1[abcd] W[aefg] p1[fehd] E2[cgbh] 
    {"abcd,efig,afhd,cgbh,ie",   4.0,   5, {11,25, 9,34,32}},  //   4.0 Ap1[abcd] W[efig] p1[afhd] E2[cgbh] delta[ie] 
    {"abcd,efig,ahfd,chbg,ie",  -4.0,   5, {11,25, 9,34,32}},  //  -4.0 Ap1[abcd] W[efig] p1[ahfd] E2[chbg] delta[ie] 
    {"abcd,abef,eghd,cgfh",      1.0,   4, {11,25,10,34}},     //   1.0 Ap1[abcd] W[abef] p2[eghd] E2[cgfh] 
    {"abcd,aefg,fbhd,cehg",      1.0,   4, {11,25,10,34}},     //   1.0 Ap1[abcd] W[aefg] p2[fbhd] E2[cehg] 
    {"abcd,aefg,fehd,cgbh",      1.0,   4, {11,25,10,34}},     //   1.0 Ap1[abcd] W[aefg] p2[fehd] E2[cgbh] 
    {"abcd,efig,afhd,cgbh,ie",  -2.0,   5, {11,25,10,34,32}},  //  -2.0 Ap1[abcd] W[efig] p2[afhd] E2[cgbh] delta[ie] 
    {"abcd,efig,ahfd,chbg,ie",   2.0,   5, {11,25,10,34,32}},  //   2.0 Ap1[abcd] W[efig] p2[ahfd] E2[chbg] delta[ie] 
    {"abcd,ebfa,eghd,cgfh",      1.0,   4, {11,26, 9,34}},     //   1.0 Ap1[abcd] W[ebfa] p1[eghd] E2[cgfh] 
    {"abcd,efga,ebhd,cfhg",      1.0,   4, {11,26, 9,34}},     //   1.0 Ap1[abcd] W[efga] p1[ebhd] E2[cfhg] 
    {"abcd,efga,eghd,cfbh",      1.0,   4, {11,26, 9,34}},     //   1.0 Ap1[abcd] W[efga] p1[eghd] E2[cfbh] 
    {"abcd,efgi,aghd,cfbh,ie",  -2.0,   5, {11,26, 9,34,32}},  //  -2.0 Ap1[abcd] W[efgi] p1[aghd] E2[cfbh] delta[ie] 
    {"abcd,efgi,ahfd,chbg,ie",   2.0,   5, {11,26, 9,34,32}},  //   2.0 Ap1[abcd] W[efgi] p1[ahfd] E2[chbg] delta[ie] 
    {"abcd,ebfa,eghd,cghf",      1.0,   4, {11,26,10,34}},     //   1.0 Ap1[abcd] W[ebfa] p2[eghd] E2[cghf] 
    {"abcd,efga,ebhd,cfgh",      1.0,   4, {11,26,10,34}},     //   1.0 Ap1[abcd] W[efga] p2[ebhd] E2[cfgh] 
    {"abcd,efga,eghd,cfbh",     -2.0,   4, {11,26,10,34}},     //  -2.0 Ap1[abcd] W[efga] p2[eghd] E2[cfbh] 
    {"abcd,efgi,aghd,cfbh,ie",   1.0,   5, {11,26,10,34,32}},  //   1.0 Ap1[abcd] W[efgi] p2[aghd] E2[cfbh] delta[ie] 
    {"abcd,efgi,ahfd,chbg,ie",  -1.0,   5, {11,26,10,34,32}},  //  -1.0 Ap1[abcd] W[efgi] p2[ahfd] E2[chbg] delta[ie] 
    {"abcd,adef,eghf,cgbh",     -2.0,   4, {11,27, 9,34}},     //  -2.0 Ap1[abcd] W[adef] p1[eghf] E2[cgbh] 
    {"abcd,edif,aghf,cgbh,ie",   4.0,   5, {11,27, 9,34,32}},  //   4.0 Ap1[abcd] W[edif] p1[aghf] E2[cgbh] delta[ie] 
    {"abcd,adef,eghf,cgbh",      1.0,   4, {11,27,10,34}},     //   1.0 Ap1[abcd] W[adef] p2[eghf] E2[cgbh] 
    {"abcd,edif,aghf,cgbh,ie",  -2.0,   5, {11,27,10,34,32}},  //  -2.0 Ap1[abcd] W[edif] p2[aghf] E2[cgbh] delta[ie] 
    {"abcd,edfa,eghf,cgbh",      4.0,   4, {11,28, 9,34}},     //   4.0 Ap1[abcd] W[edfa] p1[eghf] E2[cgbh] 
    {"abcd,edfi,aghf,cgbh,ie",  -2.0,   5, {11,28, 9,34,32}},  //  -2.0 Ap1[abcd] W[edfi] p1[aghf] E2[cgbh] delta[ie] 
    {"abcd,edfa,eghf,cgbh",     -2.0,   4, {11,28,10,34}},     //  -2.0 Ap1[abcd] W[edfa] p2[eghf] E2[cgbh] 
    {"abcd,edfi,aghf,cgbh,ie",   1.0,   5, {11,28,10,34,32}},  //   1.0 Ap1[abcd] W[edfi] p2[aghf] E2[cgbh] delta[ie] 
    {"abcd,ef,aegd,cfbg",       -1.0,   4, {12,18, 9,34}},     //  -1.0 Ap2[abcd] k[ef] p1[aegd] E2[cfbg] 
    {"abcd,ef,aged,cgbf",        1.0,   4, {12,18, 9,34}},     //   1.0 Ap2[abcd] k[ef] p1[aged] E2[cgbf] 
    {"abcd,ef,aegd,cfgb",       -1.0,   4, {12,18,10,34}},     //  -1.0 Ap2[abcd] k[ef] p2[aegd] E2[cfgb] 
    {"abcd,ef,aged,cgfb",        1.0,   4, {12,18,10,34}},     //   1.0 Ap2[abcd] k[ef] p2[aged] E2[cgfb] 
    {"abcd,ae,efgd,cfbg",        1.0,   4, {12,19, 9,34}},     //   1.0 Ap2[abcd] k[ae] p1[efgd] E2[cfbg] 
    {"abcd,ae,efgd,cfgb",        1.0,   4, {12,19,10,34}},     //   1.0 Ap2[abcd] k[ae] p2[efgd] E2[cfgb] 
    {"abcd,de,afge,cfbg",       -1.0,   4, {12,20, 9,34}},     //  -1.0 Ap2[abcd] k[de] p1[afge] E2[cfbg] 
    {"abcd,de,afge,cfgb",       -1.0,   4, {12,20,10,34}},     //  -1.0 Ap2[abcd] k[de] p2[afge] E2[cfgb] 
    {"abcd,befg,aehd,cgfh",     -1.0,   4, {12,21, 9,34}},     //  -1.0 Ap2[abcd] W[befg] p1[aehd] E2[cgfh] 
    {"abcd,befg,afhd,cehg",     -1.0,   4, {12,21, 9,34}},     //  -1.0 Ap2[abcd] W[befg] p1[afhd] E2[cehg] 
    {"abcd,befg,ahed,chfg",      1.0,   4, {12,21, 9,34}},     //   1.0 Ap2[abcd] W[befg] p1[ahed] E2[chfg] 
    {"abcd,efgh,abed,cfgh",      1.0,   4, {12,21, 9,34}},     //   1.0 Ap2[abcd] W[efgh] p1[abed] E2[cfgh] 
    {"abcd,befg,aehd,cghf",     -1.0,   4, {12,21,10,34}},     //  -1.0 Ap2[abcd] W[befg] p2[aehd] E2[cghf] 
    {"abcd,befg,afhd,cehg",      2.0,   4, {12,21,10,34}},     //   2.0 Ap2[abcd] W[befg] p2[afhd] E2[cehg] 
    {"abcd,befg,ahed,chgf",      1.0,   4, {12,21,10,34}},     //   1.0 Ap2[abcd] W[befg] p2[ahed] E2[chgf] 
    {"abcd,efgh,abed,cfgh",     -2.0,   4, {12,21,10,34}},     //  -2.0 Ap2[abcd] W[efgh] p2[abed] E2[cfgh] 
    {"abcd,bdef,aghf,cgeh",     -1.0,   4, {12,22, 9,34}},     //  -1.0 Ap2[abcd] W[bdef] p1[aghf] E2[cgeh] 
    {"abcd,edfg,abhg,cehf",     -1.0,   4, {12,22, 9,34}},     //  -1.0 Ap2[abcd] W[edfg] p1[abhg] E2[cehf] 
    {"abcd,edfg,aehg,cfbh",     -1.0,   4, {12,22, 9,34}},     //  -1.0 Ap2[abcd] W[edfg] p1[aehg] E2[cfbh] 
    {"abcd,bdef,aghf,cghe",     -1.0,   4, {12,22,10,34}},     //  -1.0 Ap2[abcd] W[bdef] p2[aghf] E2[cghe] 
    {"abcd,edfg,abhg,cehf",      2.0,   4, {12,22,10,34}},     //   2.0 Ap2[abcd] W[edfg] p2[abhg] E2[cehf] 
    {"abcd,edfg,aehg,cfhb",     -1.0,   4, {12,22,10,34}},     //  -1.0 Ap2[abcd] W[edfg] p2[aehg] E2[cfhb] 
    {"abcd,bdef,aghe,cgfh",      2.0,   4, {12,23, 9,34}},     //   2.0 Ap2[abcd] W[bdef] p1[aghe] E2[cgfh] 
    {"abcd,edfg,abhf,cegh",     -1.0,   4, {12,23, 9,34}},     //  -1.0 Ap2[abcd] W[edfg] p1[abhf] E2[cegh] 
    {"abcd,edfg,aghf,cehb",     -1.0,   4, {12,23, 9,34}},     //  -1.0 Ap2[abcd] W[edfg] p1[aghf] E2[cehb] 
    {"abcd,bdef,aghe,cgfh",     -1.0,   4, {12,23,10,34}},     //  -1.0 Ap2[abcd] W[bdef] p2[aghe] E2[cgfh] 
    {"abcd,edfg,abhf,cegh",      2.0,   4, {12,23,10,34}},     //   2.0 Ap2[abcd] W[edfg] p2[abhf] E2[cegh] 
    {"abcd,edfg,aghf,cebh",     -1.0,   4, {12,23,10,34}},     //  -1.0 Ap2[abcd] W[edfg] p2[aghf] E2[cebh] 
    {"abcd,aeif,fghd,cgbh,ie",  -1.0,   5, {12,24, 9,34,32}},  //  -1.0 Ap2[abcd] W[aeif] p1[fghd] E2[cgbh] delta[ie] 
    {"abcd,aefi,fghd,cgbh,ie",   2.0,   5, {12,24, 9,34,32}},  //   2.0 Ap2[abcd] W[aefi] p1[fghd] E2[cgbh] delta[ie] 
    {"abcd,aeif,fghd,cghb,ie",  -1.0,   5, {12,24,10,34,32}},  //  -1.0 Ap2[abcd] W[aeif] p2[fghd] E2[cghb] delta[ie] 
    {"abcd,aefi,fghd,cghb,ie",   2.0,   5, {12,24,10,34,32}},  //   2.0 Ap2[abcd] W[aefi] p2[fghd] E2[cghb] delta[ie] 
    {"abcd,abef,eghd,cgfh",      1.0,   4, {12,25, 9,34}},     //   1.0 Ap2[abcd] W[abef] p1[eghd] E2[cgfh] 
    {"abcd,aefg,fbhd,cehg",      1.0,   4, {12,25, 9,34}},     //   1.0 Ap2[abcd] W[aefg] p1[fbhd] E2[cehg] 
    {"abcd,aefg,fehd,cgbh",      1.0,   4, {12,25, 9,34}},     //   1.0 Ap2[abcd] W[aefg] p1[fehd] E2[cgbh] 
    {"abcd,efig,afhd,cgbh,ie",  -2.0,   5, {12,25, 9,34,32}},  //  -2.0 Ap2[abcd] W[efig] p1[afhd] E2[cgbh] delta[ie] 
    {"abcd,efig,ahfd,chbg,ie",   2.0,   5, {12,25, 9,34,32}},  //   2.0 Ap2[abcd] W[efig] p1[ahfd] E2[chbg] delta[ie] 
    {"abcd,abef,eghd,cghf",      1.0,   4, {12,25,10,34}},     //   1.0 Ap2[abcd] W[abef] p2[eghd] E2[cghf] 
    {"abcd,aefg,fbhd,cehg",     -2.0,   4, {12,25,10,34}},     //  -2.0 Ap2[abcd] W[aefg] p2[fbhd] E2[cehg] 
    {"abcd,aefg,fehd,cghb",      1.0,   4, {12,25,10,34}},     //   1.0 Ap2[abcd] W[aefg] p2[fehd] E2[cghb] 
    {"abcd,efig,afhd,cghb,ie",  -2.0,   5, {12,25,10,34,32}},  //  -2.0 Ap2[abcd] W[efig] p2[afhd] E2[cghb] delta[ie] 
    {"abcd,efig,ahfd,chgb,ie",   2.0,   5, {12,25,10,34,32}},  //   2.0 Ap2[abcd] W[efig] p2[ahfd] E2[chgb] delta[ie] 
    {"abcd,ebfa,eghd,cgfh",     -2.0,   4, {12,26, 9,34}},     //  -2.0 Ap2[abcd] W[ebfa] p1[eghd] E2[cgfh] 
    {"abcd,efga,ebhd,cfgh",      1.0,   4, {12,26, 9,34}},     //   1.0 Ap2[abcd] W[efga] p1[ebhd] E2[cfgh] 
    {"abcd,efga,eghd,cfhb",      1.0,   4, {12,26, 9,34}},     //   1.0 Ap2[abcd] W[efga] p1[eghd] E2[cfhb] 
    {"abcd,efgi,aghd,cfbh,ie",   1.0,   5, {12,26, 9,34,32}},  //   1.0 Ap2[abcd] W[efgi] p1[aghd] E2[cfbh] delta[ie] 
    {"abcd,efgi,ahfd,chbg,ie",  -1.0,   5, {12,26, 9,34,32}},  //  -1.0 Ap2[abcd] W[efgi] p1[ahfd] E2[chbg] delta[ie] 
    {"abcd,ebfa,eghd,cghf",     -2.0,   4, {12,26,10,34}},     //  -2.0 Ap2[abcd] W[ebfa] p2[eghd] E2[cghf] 
    {"abcd,efga,ebhd,cfhg",      1.0,   4, {12,26,10,34}},     //   1.0 Ap2[abcd] W[efga] p2[ebhd] E2[cfhg] 
    {"abcd,efga,eghd,cfhb",     -2.0,   4, {12,26,10,34}},     //  -2.0 Ap2[abcd] W[efga] p2[eghd] E2[cfhb] 
    {"abcd,efgi,aghd,cfhb,ie",   1.0,   5, {12,26,10,34,32}},  //   1.0 Ap2[abcd] W[efgi] p2[aghd] E2[cfhb] delta[ie] 
    {"abcd,efgi,ahfd,chgb,ie",  -1.0,   5, {12,26,10,34,32}},  //  -1.0 Ap2[abcd] W[efgi] p2[ahfd] E2[chgb] delta[ie] 
    {"abcd,adef,eghf,cgbh",      1.0,   4, {12,27, 9,34}},     //   1.0 Ap2[abcd] W[adef] p1[eghf] E2[cgbh] 
    {"abcd,edif,aghf,cgbh,ie",  -2.0,   5, {12,27, 9,34,32}},  //  -2.0 Ap2[abcd] W[edif] p1[aghf] E2[cgbh] delta[ie] 
    {"abcd,adef,eghf,cghb",      1.0,   4, {12,27,10,34}},     //   1.0 Ap2[abcd] W[adef] p2[eghf] E2[cghb] 
    {"abcd,edif,aghf,cghb,ie",  -2.0,   5, {12,27,10,34,32}},  //  -2.0 Ap2[abcd] W[edif] p2[aghf] E2[cghb] delta[ie] 
    {"abcd,edfa,eghf,cgbh",     -2.0,   4, {12,28, 9,34}},     //  -2.0 Ap2[abcd] W[edfa] p1[eghf] E2[cgbh] 
    {"abcd,edfi,aghf,cgbh,ie",   1.0,   5, {12,28, 9,34,32}},  //   1.0 Ap2[abcd] W[edfi] p1[aghf] E2[cgbh] delta[ie] 
    {"abcd,edfa,eghf,cgbh",      1.0,   4, {12,28,10,34}},     //   1.0 Ap2[abcd] W[edfa] p2[eghf] E2[cgbh] 
    {"abcd,edfi,aghf,cghb,ie",   1.0,   5, {12,28,10,34,32}},  //   1.0 Ap2[abcd] W[edfi] p2[aghf] E2[cghb] delta[ie] 
    {"abcd,efgh,aeid,cfgbhi",    2.0,   4, {11,21, 9,36}},     //   2.0 Ap1[abcd] W[efgh] p1[aeid] E3[cfgbhi] 
    {"abcd,efgh,aied,cfibhg",   -2.0,   4, {11,21, 9,36}},     //  -2.0 Ap1[abcd] W[efgh] p1[aied] E3[cfibhg] 
    {"abcd,efgh,aeid,cfgbhi",   -1.0,   4, {11,21,10,36}},     //  -1.0 Ap1[abcd] W[efgh] p2[aeid] E3[cfgbhi] 
    {"abcd,efgh,aied,cfibhg",    1.0,   4, {11,21,10,36}},     //   1.0 Ap1[abcd] W[efgh] p2[aied] E3[cfibhg] 
    {"abcd,edfg,ahig,cehbfi",    2.0,   4, {11,22, 9,36}},     //   2.0 Ap1[abcd] W[edfg] p1[ahig] E3[cehbfi] 
    {"abcd,edfg,ahig,cehbfi",   -1.0,   4, {11,22,10,36}},     //  -1.0 Ap1[abcd] W[edfg] p2[ahig] E3[cehbfi] 
    {"abcd,edfg,ahif,cehbgi",   -1.0,   4, {11,23, 9,36}},     //  -1.0 Ap1[abcd] W[edfg] p1[ahif] E3[cehbgi] 
    {"abcd,edfg,ahif,cehbig",   -1.0,   4, {11,23,10,36}},     //  -1.0 Ap1[abcd] W[edfg] p2[ahif] E3[cehbig] 
    {"abcd,aefg,fhid,cehbgi",   -2.0,   4, {11,25, 9,36}},     //  -2.0 Ap1[abcd] W[aefg] p1[fhid] E3[cehbgi] 
    {"abcd,aefg,fhid,cehbgi",    1.0,   4, {11,25,10,36}},     //   1.0 Ap1[abcd] W[aefg] p2[fhid] E3[cehbgi] 
    {"abcd,efga,ehid,cfhbgi",    1.0,   4, {11,26, 9,36}},     //   1.0 Ap1[abcd] W[efga] p1[ehid] E3[cfhbgi] 
    {"abcd,efga,ehid,cfhbig",    1.0,   4, {11,26,10,36}},     //   1.0 Ap1[abcd] W[efga] p2[ehid] E3[cfhbig] 
    {"abcd,efgh,aeid,cfgbhi",   -1.0,   4, {12,21, 9,36}},     //  -1.0 Ap2[abcd] W[efgh] p1[aeid] E3[cfgbhi] 
    {"abcd,efgh,aied,cfibhg",    1.0,   4, {12,21, 9,36}},     //   1.0 Ap2[abcd] W[efgh] p1[aied] E3[cfibhg] 
    {"abcd,efgh,aeid,cfgihb",   -1.0,   4, {12,21,10,36}},     //  -1.0 Ap2[abcd] W[efgh] p2[aeid] E3[cfgihb] 
    {"abcd,efgh,aied,cfighb",    1.0,   4, {12,21,10,36}},     //   1.0 Ap2[abcd] W[efgh] p2[aied] E3[cfighb] 
    {"abcd,edfg,ahig,cehbfi",   -1.0,   4, {12,22, 9,36}},     //  -1.0 Ap2[abcd] W[edfg] p1[ahig] E3[cehbfi] 
    {"abcd,edfg,ahig,cehifb",   -1.0,   4, {12,22,10,36}},     //  -1.0 Ap2[abcd] W[edfg] p2[ahig] E3[cehifb] 
    {"abcd,edfg,ahif,cehgbi",   -1.0,   4, {12,23, 9,36}},     //  -1.0 Ap2[abcd] W[edfg] p1[ahif] E3[cehgbi] 
    {"abcd,edfg,ahif,cehgib",   -1.0,   4, {12,23,10,36}},     //  -1.0 Ap2[abcd] W[edfg] p2[ahif] E3[cehgib] 
    {"abcd,aefg,fhid,cehbgi",    1.0,   4, {12,25, 9,36}},     //   1.0 Ap2[abcd] W[aefg] p1[fhid] E3[cehbgi] 
    {"abcd,aefg,fhid,cehigb",    1.0,   4, {12,25,10,36}},     //   1.0 Ap2[abcd] W[aefg] p2[fhid] E3[cehigb] 
    {"abcd,efga,ehid,cfhgbi",    1.0,   4, {12,26, 9,36}},     //   1.0 Ap2[abcd] W[efga] p1[ehid] E3[cfhgbi] 
    {"abcd,efga,ehid,cfhibg",    1.0,   4, {12,26,10,36}},     //   1.0 Ap2[abcd] W[efga] p2[ehid] E3[cfhibg] 
	};

  FEqInfo bVec[8] = {
    {"IPQA,RQSP,ARIS",   2.0,   3, {13,34,29}}, //  2.0 b1[IPQA] E2[RQSP]           V1[ARIS]
    {"IPQA,SQ,APIS"  ,   2.0,   3, {13,33,29}}, //  2.0 b1[IPQA] E1[SQ]   delta[PR] V1[ARIS]
    {"IPQA,RQSP,RAIS",  -1.0,   3, {13,34,30}}, // -1.0 b1[IPQA] E2[RQSP]           V2[RAIS]
    {"IPQA,SQ,PAIS"  ,  -1.0,   3, {13,33,30}}, // -1.0 b1[IPQA] E1[SQ]   delta[PR] V2[RAIS]
    {"IPQA,RQSP,ARIS",  -1.0,   3, {14,34,29}}, // -1.0 b2[IPQA] E2[RQSP]           V1[ARIS]
    {"IPQA,SQ,APIS"  ,  -1.0,   3, {14,33,29}}, // -1.0 b2[IPQA] E1[SQ]   delta[PR] V1[ARIS]
    {"IPQA,RQPS,RAIS",  -1.0,   3, {14,34,30}}, // -1.0 b2[IPQA] E2[RQPS]           V2[RAIS]
    {"IPQA,SQ,PAIS"  ,   2.0,   3, {14,33,30}}, //  2.0 b2[IPQA] E1[SQ]   delta[PR] V2[RAIS]
  };

  int f(int i) {
    return 2*i;
  }
  FDomainDecl DomainDecls[1] = {
    {"A", "a", f}
  };

  static void GetMethodInfo(FMethodInfo &Out) {
    Out = FMethodInfo();
    Out.pName = "MRLCC_CAAV";
    Out.perturberClass = "CAAV";
    Out.pSpinClass = "restricted";
    Out.pTensorDecls = &TensorDecls[0];
    Out.nTensorDecls = 37;
    Out.pDomainDecls = &DomainDecls[0];
    Out.nDomainDecls = 1;
    Out.EqsRes = FEqSet(&EqsRes[0], 216, "MRLCC_CAAV/Res");
    Out.bVec = FEqSet(&bVec[0], 8, "MRLCC_CAAV/bVec");
  };
};
