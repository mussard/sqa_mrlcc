namespace MRLCC_CCAV {

  FTensorDecl TensorDecls[27] = {
    /*  0*/{"R"      , "ccae"     , "", USAGE_Residual    },
    /*  1*/{"t"      , "ccae"     , "", USAGE_Amplitude   },
    /*  2*/{"T"      , "ccae"     , "", USAGE_Amplitude   },
    /*  3*/{"b"      , "ccae"     , "", USAGE_Amplitude   },
    /*  4*/{"B"      , "ccae"     , "", USAGE_Amplitude   },
    /*  5*/{"p"      , "ccae"     , "", USAGE_Amplitude   },
    /*  6*/{"Ap"     , "ccae"     , "", USAGE_Amplitude   },
    /*  7*/{"P"      , "ccae"     , "", USAGE_Amplitude   },
    /*  8*/{"AP"     , "ccae"     , "", USAGE_Amplitude   },
    /*  9*/{"k"      , "aa"       , "", USAGE_Hamiltonian },
    /* 10*/{"k"      , "cc"       , "", USAGE_Hamiltonian },
    /* 11*/{"k"      , "ee"       , "", USAGE_Hamiltonian },
    /* 12*/{"W"      , "aaaa"     , "", USAGE_Hamiltonian },
    /* 13*/{"W"      , "aeae"     , "", USAGE_Hamiltonian },
    /* 14*/{"W"      , "aeea"     , "", USAGE_Hamiltonian },
    /* 15*/{"W"      , "cccc"     , "", USAGE_Hamiltonian },
    /* 16*/{"W"      , "caca"     , "", USAGE_Hamiltonian },
    /* 17*/{"W"      , "caac"     , "", USAGE_Hamiltonian },
    /* 18*/{"W"      , "cece"     , "", USAGE_Hamiltonian },
    /* 19*/{"W"      , "ceec"     , "", USAGE_Hamiltonian },
    /* 20*/{"W"      , "ccae"     , "", USAGE_Hamiltonian },
    /* 21*/{"delta"  , "aa"       , "", USAGE_Density     },
    /* 22*/{"delta"  , "cc"       , "", USAGE_Density     },
    /* 23*/{"E1"     , "aa"       , "", USAGE_Density     },
    /* 24*/{"E2"     , "aaaa"     , "", USAGE_Density     },
    /* 25*/{"S1"     , "AA"       , "", USAGE_Density     },
    /* 26*/{"S2"     , "aa"       , "", USAGE_Density     },
  };

	FEqInfo EqsRes[134] = {
    {"abcd,ce,abed",           4.0,   3, { 6, 9, 5}},        //   4.0 Ap[abcd] k[ce] p[abed] 
    {"abcd,ce,baed",          -2.0,   3, { 6, 9, 5}},        //  -2.0 Ap[abcd] k[ce] p[baed] 
    {"abcd,ae,becd",           2.0,   3, { 6,10, 5}},        //   2.0 Ap[abcd] k[ae] p[becd] 
    {"abcd,ae,ebcd",          -4.0,   3, { 6,10, 5}},        //  -4.0 Ap[abcd] k[ae] p[ebcd] 
    {"abcd,be,aecd",          -4.0,   3, { 6,10, 5}},        //  -4.0 Ap[abcd] k[be] p[aecd] 
    {"abcd,be,eacd",           2.0,   3, { 6,10, 5}},        //   2.0 Ap[abcd] k[be] p[eacd] 
    {"abcd,de,abce",           4.0,   3, { 6,11, 5}},        //   4.0 Ap[abcd] k[de] p[abce] 
    {"abcd,de,bace",          -2.0,   3, { 6,11, 5}},        //  -2.0 Ap[abcd] k[de] p[bace] 
    {"abcd,cdef,abef",         4.0,   3, { 6,13, 5}},        //   4.0 Ap[abcd] W[cdef] p[abef] 
    {"abcd,cdef,baef",        -2.0,   3, { 6,13, 5}},        //  -2.0 Ap[abcd] W[cdef] p[baef] 
    {"abcd,cdef,abfe",        -2.0,   3, { 6,14, 5}},        //  -2.0 Ap[abcd] W[cdef] p[abfe] 
    {"abcd,cdef,bafe",         4.0,   3, { 6,14, 5}},        //   4.0 Ap[abcd] W[cdef] p[bafe] 
    {"abcd,abef,efcd",         4.0,   3, { 6,15, 5}},        //   4.0 Ap[abcd] W[abef] p[efcd] 
    {"abcd,abef,fecd",        -2.0,   3, { 6,15, 5}},        //  -2.0 Ap[abcd] W[abef] p[fecd] 
    {"abcd,aegf,bfcd,ge",     -2.0,   4, { 6,15, 5,22}},     //  -2.0 Ap[abcd] W[aegf] p[bfcd] delta[ge] 
    {"abcd,aegf,fbcd,ge",      4.0,   4, { 6,15, 5,22}},     //   4.0 Ap[abcd] W[aegf] p[fbcd] delta[ge] 
    {"abcd,aefg,bfcd,ge",      4.0,   4, { 6,15, 5,22}},     //   4.0 Ap[abcd] W[aefg] p[bfcd] delta[ge] 
    {"abcd,aefg,fbcd,ge",     -8.0,   4, { 6,15, 5,22}},     //  -8.0 Ap[abcd] W[aefg] p[fbcd] delta[ge] 
    {"abcd,begf,afcd,ge",      4.0,   4, { 6,15, 5,22}},     //   4.0 Ap[abcd] W[begf] p[afcd] delta[ge] 
    {"abcd,begf,facd,ge",     -2.0,   4, { 6,15, 5,22}},     //  -2.0 Ap[abcd] W[begf] p[facd] delta[ge] 
    {"abcd,befg,afcd,ge",     -8.0,   4, { 6,15, 5,22}},     //  -8.0 Ap[abcd] W[befg] p[afcd] delta[ge] 
    {"abcd,befg,facd,ge",      4.0,   4, { 6,15, 5,22}},     //   4.0 Ap[abcd] W[befg] p[facd] delta[ge] 
    {"abcd,acef,befd",         2.0,   3, { 6,16, 5}},        //   2.0 Ap[abcd] W[acef] p[befd] 
    {"abcd,acef,ebfd",        -4.0,   3, { 6,16, 5}},        //  -4.0 Ap[abcd] W[acef] p[ebfd] 
    {"abcd,bcef,aefd",        -4.0,   3, { 6,16, 5}},        //  -4.0 Ap[abcd] W[bcef] p[aefd] 
    {"abcd,bcef,eafd",         2.0,   3, { 6,16, 5}},        //   2.0 Ap[abcd] W[bcef] p[eafd] 
    {"abcd,ecgf,abfd,ge",      8.0,   4, { 6,16, 5,22}},     //   8.0 Ap[abcd] W[ecgf] p[abfd] delta[ge] 
    {"abcd,ecgf,bafd,ge",     -4.0,   4, { 6,16, 5,22}},     //  -4.0 Ap[abcd] W[ecgf] p[bafd] delta[ge] 
    {"abcd,ecfa,befd",        -4.0,   3, { 6,17, 5}},        //  -4.0 Ap[abcd] W[ecfa] p[befd] 
    {"abcd,ecfa,ebfd",         8.0,   3, { 6,17, 5}},        //   8.0 Ap[abcd] W[ecfa] p[ebfd] 
    {"abcd,ecfb,aefd",         2.0,   3, { 6,17, 5}},        //   2.0 Ap[abcd] W[ecfb] p[aefd] 
    {"abcd,ecfb,eafd",        -4.0,   3, { 6,17, 5}},        //  -4.0 Ap[abcd] W[ecfb] p[eafd] 
    {"abcd,ecfg,abfd,ge",     -4.0,   4, { 6,17, 5,22}},     //  -4.0 Ap[abcd] W[ecfg] p[abfd] delta[ge] 
    {"abcd,ecfg,bafd,ge",      2.0,   4, { 6,17, 5,22}},     //   2.0 Ap[abcd] W[ecfg] p[bafd] delta[ge] 
    {"abcd,adef,becf",         2.0,   3, { 6,18, 5}},        //   2.0 Ap[abcd] W[adef] p[becf] 
    {"abcd,adef,ebcf",        -4.0,   3, { 6,18, 5}},        //  -4.0 Ap[abcd] W[adef] p[ebcf] 
    {"abcd,bdef,aecf",        -4.0,   3, { 6,18, 5}},        //  -4.0 Ap[abcd] W[bdef] p[aecf] 
    {"abcd,bdef,eacf",         2.0,   3, { 6,18, 5}},        //   2.0 Ap[abcd] W[bdef] p[eacf] 
    {"abcd,edgf,abcf,ge",      8.0,   4, { 6,18, 5,22}},     //   8.0 Ap[abcd] W[edgf] p[abcf] delta[ge] 
    {"abcd,edgf,bacf,ge",     -4.0,   4, { 6,18, 5,22}},     //  -4.0 Ap[abcd] W[edgf] p[bacf] delta[ge] 
    {"abcd,edfa,becf",        -4.0,   3, { 6,19, 5}},        //  -4.0 Ap[abcd] W[edfa] p[becf] 
    {"abcd,edfa,ebcf",         2.0,   3, { 6,19, 5}},        //   2.0 Ap[abcd] W[edfa] p[ebcf] 
    {"abcd,edfb,aecf",         8.0,   3, { 6,19, 5}},        //   8.0 Ap[abcd] W[edfb] p[aecf] 
    {"abcd,edfb,eacf",        -4.0,   3, { 6,19, 5}},        //  -4.0 Ap[abcd] W[edfb] p[eacf] 
    {"abcd,edfg,abcf,ge",     -4.0,   4, { 6,19, 5,22}},     //  -4.0 Ap[abcd] W[edfg] p[abcf] delta[ge] 
    {"abcd,edfg,bacf,ge",      2.0,   4, { 6,19, 5,22}},     //   2.0 Ap[abcd] W[edfg] p[bacf] delta[ge] 
    {"abcd,ef,abed,fc",       -2.0,   4, { 6, 9, 5,23}},     //  -2.0 Ap[abcd] k[ef] p[abed] E1[fc] 
    {"abcd,ef,baed,fc",        1.0,   4, { 6, 9, 5,23}},     //   1.0 Ap[abcd] k[ef] p[baed] E1[fc] 
    {"abcd,ae,befd,fc",       -1.0,   4, { 6,10, 5,23}},     //  -1.0 Ap[abcd] k[ae] p[befd] E1[fc] 
    {"abcd,ae,ebfd,fc",        2.0,   4, { 6,10, 5,23}},     //   2.0 Ap[abcd] k[ae] p[ebfd] E1[fc] 
    {"abcd,be,aefd,fc",        2.0,   4, { 6,10, 5,23}},     //   2.0 Ap[abcd] k[be] p[aefd] E1[fc] 
    {"abcd,be,eafd,fc",       -1.0,   4, { 6,10, 5,23}},     //  -1.0 Ap[abcd] k[be] p[eafd] E1[fc] 
    {"abcd,de,abfe,fc",       -2.0,   4, { 6,11, 5,23}},     //  -2.0 Ap[abcd] k[de] p[abfe] E1[fc] 
    {"abcd,de,bafe,fc",        1.0,   4, { 6,11, 5,23}},     //   1.0 Ap[abcd] k[de] p[bafe] E1[fc] 
    {"abcd,cefg,abed,gf",     -2.0,   4, { 6,12, 5,23}},     //  -2.0 Ap[abcd] W[cefg] p[abed] E1[gf] 
    {"abcd,cefg,abfd,eg",      4.0,   4, { 6,12, 5,23}},     //   4.0 Ap[abcd] W[cefg] p[abfd] E1[eg] 
    {"abcd,cefg,baed,gf",      1.0,   4, { 6,12, 5,23}},     //   1.0 Ap[abcd] W[cefg] p[baed] E1[gf] 
    {"abcd,cefg,bafd,eg",     -2.0,   4, { 6,12, 5,23}},     //  -2.0 Ap[abcd] W[cefg] p[bafd] E1[eg] 
    {"abcd,cdef,abgf,ge",     -2.0,   4, { 6,13, 5,23}},     //  -2.0 Ap[abcd] W[cdef] p[abgf] E1[ge] 
    {"abcd,cdef,bagf,ge",      1.0,   4, { 6,13, 5,23}},     //   1.0 Ap[abcd] W[cdef] p[bagf] E1[ge] 
    {"abcd,edfg,abcg,ef",      4.0,   4, { 6,13, 5,23}},     //   4.0 Ap[abcd] W[edfg] p[abcg] E1[ef] 
    {"abcd,edfg,abeg,fc",     -2.0,   4, { 6,13, 5,23}},     //  -2.0 Ap[abcd] W[edfg] p[abeg] E1[fc] 
    {"abcd,edfg,bacg,ef",     -2.0,   4, { 6,13, 5,23}},     //  -2.0 Ap[abcd] W[edfg] p[bacg] E1[ef] 
    {"abcd,edfg,baeg,fc",      1.0,   4, { 6,13, 5,23}},     //   1.0 Ap[abcd] W[edfg] p[baeg] E1[fc] 
    {"abcd,cdef,abge,gf",      1.0,   4, { 6,14, 5,23}},     //   1.0 Ap[abcd] W[cdef] p[abge] E1[gf] 
    {"abcd,cdef,bage,gf",     -2.0,   4, { 6,14, 5,23}},     //  -2.0 Ap[abcd] W[cdef] p[bage] E1[gf] 
    {"abcd,edfg,abcf,eg",     -2.0,   4, { 6,14, 5,23}},     //  -2.0 Ap[abcd] W[edfg] p[abcf] E1[eg] 
    {"abcd,edfg,abgf,ec",      1.0,   4, { 6,14, 5,23}},     //   1.0 Ap[abcd] W[edfg] p[abgf] E1[ec] 
    {"abcd,edfg,bacf,eg",      1.0,   4, { 6,14, 5,23}},     //   1.0 Ap[abcd] W[edfg] p[bacf] E1[eg] 
    {"abcd,edfg,bagf,ec",     -2.0,   4, { 6,14, 5,23}},     //  -2.0 Ap[abcd] W[edfg] p[bagf] E1[ec] 
    {"abcd,abef,efgd,gc",     -2.0,   4, { 6,15, 5,23}},     //  -2.0 Ap[abcd] W[abef] p[efgd] E1[gc] 
    {"abcd,abef,fegd,gc",      1.0,   4, { 6,15, 5,23}},     //   1.0 Ap[abcd] W[abef] p[fegd] E1[gc] 
    {"abcd,aehf,bfgd,gc,he",   1.0,   5, { 6,15, 5,23,22}},  //   1.0 Ap[abcd] W[aehf] p[bfgd] E1[gc] delta[he] 
    {"abcd,aehf,fbgd,gc,he",  -2.0,   5, { 6,15, 5,23,22}},  //  -2.0 Ap[abcd] W[aehf] p[fbgd] E1[gc] delta[he] 
    {"abcd,aefh,bfgd,gc,he",  -2.0,   5, { 6,15, 5,23,22}},  //  -2.0 Ap[abcd] W[aefh] p[bfgd] E1[gc] delta[he] 
    {"abcd,aefh,fbgd,gc,he",   4.0,   5, { 6,15, 5,23,22}},  //   4.0 Ap[abcd] W[aefh] p[fbgd] E1[gc] delta[he] 
    {"abcd,behf,afgd,gc,he",  -2.0,   5, { 6,15, 5,23,22}},  //  -2.0 Ap[abcd] W[behf] p[afgd] E1[gc] delta[he] 
    {"abcd,behf,fagd,gc,he",   1.0,   5, { 6,15, 5,23,22}},  //   1.0 Ap[abcd] W[behf] p[fagd] E1[gc] delta[he] 
    {"abcd,befh,afgd,gc,he",   4.0,   5, { 6,15, 5,23,22}},  //   4.0 Ap[abcd] W[befh] p[afgd] E1[gc] delta[he] 
    {"abcd,befh,fagd,gc,he",  -2.0,   5, { 6,15, 5,23,22}},  //  -2.0 Ap[abcd] W[befh] p[fagd] E1[gc] delta[he] 
    {"abcd,acef,begd,gf",     -1.0,   4, { 6,16, 5,23}},     //  -1.0 Ap[abcd] W[acef] p[begd] E1[gf] 
    {"abcd,acef,ebgd,gf",      2.0,   4, { 6,16, 5,23}},     //   2.0 Ap[abcd] W[acef] p[ebgd] E1[gf] 
    {"abcd,aefg,bfcd,eg",      2.0,   4, { 6,16, 5,23}},     //   2.0 Ap[abcd] W[aefg] p[bfcd] E1[eg] 
    {"abcd,aefg,bfed,gc",     -1.0,   4, { 6,16, 5,23}},     //  -1.0 Ap[abcd] W[aefg] p[bfed] E1[gc] 
    {"abcd,aefg,fbcd,eg",     -4.0,   4, { 6,16, 5,23}},     //  -4.0 Ap[abcd] W[aefg] p[fbcd] E1[eg] 
    {"abcd,aefg,fbed,gc",      2.0,   4, { 6,16, 5,23}},     //   2.0 Ap[abcd] W[aefg] p[fbed] E1[gc] 
    {"abcd,bcef,aegd,gf",      2.0,   4, { 6,16, 5,23}},     //   2.0 Ap[abcd] W[bcef] p[aegd] E1[gf] 
    {"abcd,bcef,eagd,gf",     -1.0,   4, { 6,16, 5,23}},     //  -1.0 Ap[abcd] W[bcef] p[eagd] E1[gf] 
    {"abcd,befg,afcd,eg",     -4.0,   4, { 6,16, 5,23}},     //  -4.0 Ap[abcd] W[befg] p[afcd] E1[eg] 
    {"abcd,befg,afed,gc",      2.0,   4, { 6,16, 5,23}},     //   2.0 Ap[abcd] W[befg] p[afed] E1[gc] 
    {"abcd,befg,facd,eg",      2.0,   4, { 6,16, 5,23}},     //   2.0 Ap[abcd] W[befg] p[facd] E1[eg] 
    {"abcd,befg,faed,gc",     -1.0,   4, { 6,16, 5,23}},     //  -1.0 Ap[abcd] W[befg] p[faed] E1[gc] 
    {"abcd,efhg,abfd,gc,he",  -4.0,   5, { 6,16, 5,23,22}},  //  -4.0 Ap[abcd] W[efhg] p[abfd] E1[gc] delta[he] 
    {"abcd,efhg,bafd,gc,he",   2.0,   5, { 6,16, 5,23,22}},  //   2.0 Ap[abcd] W[efhg] p[bafd] E1[gc] delta[he] 
    {"abcd,ecfa,begd,gf",      2.0,   4, { 6,17, 5,23}},     //   2.0 Ap[abcd] W[ecfa] p[begd] E1[gf] 
    {"abcd,ecfa,ebgd,gf",     -4.0,   4, { 6,17, 5,23}},     //  -4.0 Ap[abcd] W[ecfa] p[ebgd] E1[gf] 
    {"abcd,ecfb,aegd,gf",     -1.0,   4, { 6,17, 5,23}},     //  -1.0 Ap[abcd] W[ecfb] p[aegd] E1[gf] 
    {"abcd,ecfb,eagd,gf",      2.0,   4, { 6,17, 5,23}},     //   2.0 Ap[abcd] W[ecfb] p[eagd] E1[gf] 
    {"abcd,efga,becd,fg",     -1.0,   4, { 6,17, 5,23}},     //  -1.0 Ap[abcd] W[efga] p[becd] E1[fg] 
    {"abcd,efga,begd,fc",      2.0,   4, { 6,17, 5,23}},     //   2.0 Ap[abcd] W[efga] p[begd] E1[fc] 
    {"abcd,efga,ebcd,fg",      2.0,   4, { 6,17, 5,23}},     //   2.0 Ap[abcd] W[efga] p[ebcd] E1[fg] 
    {"abcd,efga,ebgd,fc",     -4.0,   4, { 6,17, 5,23}},     //  -4.0 Ap[abcd] W[efga] p[ebgd] E1[fc] 
    {"abcd,efgb,aecd,fg",      2.0,   4, { 6,17, 5,23}},     //   2.0 Ap[abcd] W[efgb] p[aecd] E1[fg] 
    {"abcd,efgb,aegd,fc",     -1.0,   4, { 6,17, 5,23}},     //  -1.0 Ap[abcd] W[efgb] p[aegd] E1[fc] 
    {"abcd,efgb,eacd,fg",     -1.0,   4, { 6,17, 5,23}},     //  -1.0 Ap[abcd] W[efgb] p[eacd] E1[fg] 
    {"abcd,efgb,eagd,fc",      2.0,   4, { 6,17, 5,23}},     //   2.0 Ap[abcd] W[efgb] p[eagd] E1[fc] 
    {"abcd,efgh,abgd,fc,he",   2.0,   5, { 6,17, 5,23,22}},  //   2.0 Ap[abcd] W[efgh] p[abgd] E1[fc] delta[he] 
    {"abcd,efgh,bagd,fc,he",  -1.0,   5, { 6,17, 5,23,22}},  //  -1.0 Ap[abcd] W[efgh] p[bagd] E1[fc] delta[he] 
    {"abcd,adef,begf,gc",     -1.0,   4, { 6,18, 5,23}},     //  -1.0 Ap[abcd] W[adef] p[begf] E1[gc] 
    {"abcd,adef,ebgf,gc",      2.0,   4, { 6,18, 5,23}},     //   2.0 Ap[abcd] W[adef] p[ebgf] E1[gc] 
    {"abcd,bdef,aegf,gc",      2.0,   4, { 6,18, 5,23}},     //   2.0 Ap[abcd] W[bdef] p[aegf] E1[gc] 
    {"abcd,bdef,eagf,gc",     -1.0,   4, { 6,18, 5,23}},     //  -1.0 Ap[abcd] W[bdef] p[eagf] E1[gc] 
    {"abcd,edhf,abgf,gc,he",  -4.0,   5, { 6,18, 5,23,22}},  //  -4.0 Ap[abcd] W[edhf] p[abgf] E1[gc] delta[he] 
    {"abcd,edhf,bagf,gc,he",   2.0,   5, { 6,18, 5,23,22}},  //   2.0 Ap[abcd] W[edhf] p[bagf] E1[gc] delta[he] 
    {"abcd,edfa,begf,gc",      2.0,   4, { 6,19, 5,23}},     //   2.0 Ap[abcd] W[edfa] p[begf] E1[gc] 
    {"abcd,edfa,ebgf,gc",     -1.0,   4, { 6,19, 5,23}},     //  -1.0 Ap[abcd] W[edfa] p[ebgf] E1[gc] 
    {"abcd,edfb,aegf,gc",     -4.0,   4, { 6,19, 5,23}},     //  -4.0 Ap[abcd] W[edfb] p[aegf] E1[gc] 
    {"abcd,edfb,eagf,gc",      2.0,   4, { 6,19, 5,23}},     //   2.0 Ap[abcd] W[edfb] p[eagf] E1[gc] 
    {"abcd,edfh,abgf,gc,he",   2.0,   5, { 6,19, 5,23,22}},  //   2.0 Ap[abcd] W[edfh] p[abgf] E1[gc] delta[he] 
    {"abcd,edfh,bagf,gc,he",  -1.0,   5, { 6,19, 5,23,22}},  //  -1.0 Ap[abcd] W[edfh] p[bagf] E1[gc] delta[he] 
    {"abcd,efgh,abed,fghc",   -2.0,   4, { 6,12, 5,24}},     //  -2.0 Ap[abcd] W[efgh] p[abed] E2[fghc] 
    {"abcd,efgh,baed,fghc",    1.0,   4, { 6,12, 5,24}},     //   1.0 Ap[abcd] W[efgh] p[baed] E2[fghc] 
    {"abcd,edfg,abhg,ehfc",   -2.0,   4, { 6,13, 5,24}},     //  -2.0 Ap[abcd] W[edfg] p[abhg] E2[ehfc] 
    {"abcd,edfg,bahg,ehfc",    1.0,   4, { 6,13, 5,24}},     //   1.0 Ap[abcd] W[edfg] p[bahg] E2[ehfc] 
    {"abcd,edfg,abhf,ehgc",    1.0,   4, { 6,14, 5,24}},     //   1.0 Ap[abcd] W[edfg] p[abhf] E2[ehgc] 
    {"abcd,edfg,bahf,ehcg",    1.0,   4, { 6,14, 5,24}},     //   1.0 Ap[abcd] W[edfg] p[bahf] E2[ehcg] 
    {"abcd,aefg,bfhd,ehgc",   -1.0,   4, { 6,16, 5,24}},     //  -1.0 Ap[abcd] W[aefg] p[bfhd] E2[ehgc] 
    {"abcd,aefg,fbhd,ehgc",    2.0,   4, { 6,16, 5,24}},     //   2.0 Ap[abcd] W[aefg] p[fbhd] E2[ehgc] 
    {"abcd,befg,afhd,ehgc",    2.0,   4, { 6,16, 5,24}},     //   2.0 Ap[abcd] W[befg] p[afhd] E2[ehgc] 
    {"abcd,befg,fahd,ehgc",   -1.0,   4, { 6,16, 5,24}},     //  -1.0 Ap[abcd] W[befg] p[fahd] E2[ehgc] 
    {"abcd,efga,behd,fhcg",   -1.0,   4, { 6,17, 5,24}},     //  -1.0 Ap[abcd] W[efga] p[behd] E2[fhcg] 
    {"abcd,efga,ebhd,fhcg",    2.0,   4, { 6,17, 5,24}},     //   2.0 Ap[abcd] W[efga] p[ebhd] E2[fhcg] 
    {"abcd,efgb,aehd,fhgc",   -1.0,   4, { 6,17, 5,24}},     //  -1.0 Ap[abcd] W[efgb] p[aehd] E2[fhgc] 
    {"abcd,efgb,eahd,fhcg",   -1.0,   4, { 6,17, 5,24}},     //  -1.0 Ap[abcd] W[efgb] p[eahd] E2[fhcg] 
	};

  FEqInfo bVec[4] = {
    {"IJPA,RP,IJRA",   4.0,   3, { 3,21,20}}, //  4.0 b[IJPA] delta[RP] V[IJRA]
    {"IJPA,RP,IJRA",  -2.0,   3, { 3,23,20}}, // -2.0 b[IJPA] E1[RP]    V[IJRA]
    {"IJPA,RP,JIRA",  -2.0,   3, { 3,21,20}}, // -2.0 b[IJPA] delta[RP] V[JIRA]
    {"IJPA,RP,JIRA",   1.0,   3, { 3,23,20}}, //  1.0 b[IJPA] E1[RP]    V[JIRA]
  };

  int f(int i) {
    return 2*i;
  }
  FDomainDecl DomainDecls[1] = {
    {"A", "a", f}
  };

  static void GetMethodInfo(FMethodInfo &Out) {
    Out = FMethodInfo();
    Out.pName = "MRLCC_CCAV";
    Out.perturberClass = "CCAV";
    Out.pSpinClass = "restricted";
    Out.pTensorDecls = &TensorDecls[0];
    Out.nTensorDecls = 27;
    Out.pDomainDecls = &DomainDecls[0];
    Out.nDomainDecls = 1;
    Out.EqsRes = FEqSet(&EqsRes[0], 134, "MRLCC_CCAV/Res");
    Out.bVec = FEqSet(&bVec[0], 4, "MRLCC_CCAV/bVec");
  };
};
