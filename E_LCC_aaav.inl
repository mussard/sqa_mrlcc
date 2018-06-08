namespace MRLCC_AAAV {

  FTensorDecl TensorDecls[27] = {
    /*  0*/{"R"      , "eaaa"     , "", USAGE_Residual    },
    /*  1*/{"t"      , "eaaa"     , "", USAGE_Amplitude   },
    /*  2*/{"T"      , "eaaa"     , "", USAGE_Amplitude   },
    /*  3*/{"b"      , "eaaa"     , "", USAGE_Amplitude   },
    /*  4*/{"B"      , "eaaa"     , "", USAGE_Amplitude   },
    /*  5*/{"p"      , "eaaa"     , "", USAGE_Amplitude   },
    /*  6*/{"Ap"     , "eaaa"     , "", USAGE_Amplitude   },
    /*  7*/{"P"      , "eaaa"     , "", USAGE_Amplitude   },
    /*  8*/{"AP"     , "eaaa"     , "", USAGE_Amplitude   },
    /*  9*/{"f"      , "ea"       , "", USAGE_Hamiltonian },
    /* 10*/{"k"      , "aa"       , "", USAGE_Hamiltonian },
    /* 11*/{"k"      , "ee"       , "", USAGE_Hamiltonian },
    /* 12*/{"W"      , "aaaa"     , "", USAGE_Hamiltonian },
    /* 13*/{"W"      , "aeae"     , "", USAGE_Hamiltonian },
    /* 14*/{"W"      , "aeea"     , "", USAGE_Hamiltonian },
    /* 15*/{"W"      , "caca"     , "", USAGE_Hamiltonian },
    /* 16*/{"W"      , "caac"     , "", USAGE_Hamiltonian },
    /* 17*/{"W"      , "cece"     , "", USAGE_Hamiltonian },
    /* 18*/{"W"      , "ceec"     , "", USAGE_Hamiltonian },
    /* 19*/{"W"      , "eaaa"     , "", USAGE_Hamiltonian },
    /* 20*/{"delta"  , "aa"       , "", USAGE_Density     },
    /* 21*/{"delta"  , "cc"       , "", USAGE_Density     },
    /* 22*/{"E1"     , "aa"       , "", USAGE_Density     },
    /* 23*/{"E2"     , "aaaa"     , "", USAGE_Density     },
    /* 24*/{"S1"     , "aaaaaa"   , "", USAGE_Density     },
    /* 25*/{"E3"     , "aaaaaa"   , "", USAGE_Density     },
    /* 26*/{"E4"     , "aaaaaaaa" , "", USAGE_Density     },
  };

	FEqInfo EqsRes[86] = {
    {"abgc,de,afhd,ce,gb,hf",       -1.0,   6, { 6,10, 5,22,20,20}},     //  -1.0 Ap[abgc] k[de] p[afhd] E1[ce] delta[gb] delta[hf] 
    {"abgc,ad,dehf,cf,gb,he",        1.0,   6, { 6,11, 5,22,20,20}},     //   1.0 Ap[abgc] k[ad] p[dehf] E1[cf] delta[gb] delta[he] 
    {"abhc,deif,agje,cf,hb,id,jg",  -2.0,   7, { 6,15, 5,22,20,21,20}},  //  -2.0 Ap[abhc] W[deif] p[agje] E1[cf] delta[hb] delta[id] delta[jg] 
    {"abhc,defi,agje,cf,hb,id,jg",   1.0,   7, { 6,16, 5,22,20,21,20}},  //   1.0 Ap[abhc] W[defi] p[agje] E1[cf] delta[hb] delta[id] delta[jg] 
    {"abhc,daie,efjg,cg,hb,id,jf",   2.0,   7, { 6,17, 5,22,20,21,20}},  //   2.0 Ap[abhc] W[daie] p[efjg] E1[cg] delta[hb] delta[id] delta[jf] 
    {"abhc,daei,efjg,cg,hb,id,jf",  -1.0,   7, { 6,18, 5,22,20,21,20}},  //  -1.0 Ap[abhc] W[daei] p[efjg] E1[cg] delta[hb] delta[id] delta[jf] 
    {"abhc,de,adfg,cefg,hb",         1.0,   5, { 6,10, 5,23,20}},        //   1.0 Ap[abhc] k[de] p[adfg] E2[cefg] delta[hb] 
    {"abhc,de,afdg,cfeg,hb",        -1.0,   5, { 6,10, 5,23,20}},        //  -1.0 Ap[abhc] k[de] p[afdg] E2[cfeg] delta[hb] 
    {"abhc,de,afgd,cfge,hb",        -1.0,   5, { 6,10, 5,23,20}},        //  -1.0 Ap[abhc] k[de] p[afgd] E2[cfge] delta[hb] 
    {"abcd,be,aefg,cdfg",            1.0,   4, { 6,10, 5,23}},           //   1.0 Ap[abcd] k[be] p[aefg] E2[cdfg] 
    {"abcd,ef,abeg,cdfg",           -1.0,   4, { 6,10, 5,23}},           //  -1.0 Ap[abcd] k[ef] p[abeg] E2[cdfg] 
    {"abcd,ef,abge,cdgf",           -1.0,   4, { 6,10, 5,23}},           //  -1.0 Ap[abcd] k[ef] p[abge] E2[cdgf] 
    {"abcd,ef,aghe,cdfb,hg",        -1.0,   5, { 6,10, 5,23,20}},        //  -1.0 Ap[abcd] k[ef] p[aghe] E2[cdfb] delta[hg] 
    {"abhc,ad,defg,cefg,hb",         1.0,   5, { 6,11, 5,23,20}},        //   1.0 Ap[abhc] k[ad] p[defg] E2[cefg] delta[hb] 
    {"abcd,ae,ebfg,cdfg",            1.0,   4, { 6,11, 5,23}},           //   1.0 Ap[abcd] k[ae] p[ebfg] E2[cdfg] 
    {"abcd,ae,efhg,cdgb,hf",         1.0,   5, { 6,11, 5,23,20}},        //   1.0 Ap[abcd] k[ae] p[efhg] E2[cdgb] delta[hf] 
    {"abic,defg,ahde,chfg,ib",      -1.0,   5, { 6,12, 5,23,20}},        //  -1.0 Ap[abic] W[defg] p[ahde] E2[chfg] delta[ib] 
    {"abic,defg,ahjd,cefg,ib,jh",   -1.0,   6, { 6,12, 5,23,20,20}},     //  -1.0 Ap[abic] W[defg] p[ahjd] E2[cefg] delta[ib] delta[jh] 
    {"abcd,befg,ahie,cdgf,ih",      -1.0,   5, { 6,12, 5,23,20}},        //  -1.0 Ap[abcd] W[befg] p[ahie] E2[cdgf] delta[ih] 
    {"abcd,efgh,abef,cdgh",         -1.0,   4, { 6,12, 5,23}},           //  -1.0 Ap[abcd] W[efgh] p[abef] E2[cdgh] 
    {"abic,daef,fdgh,cegh,ib",       1.0,   5, { 6,13, 5,23,20}},        //   1.0 Ap[abic] W[daef] p[fdgh] E2[cegh] delta[ib] 
    {"abic,daef,fgjh,cdhe,ib,jg",    1.0,   6, { 6,13, 5,23,20,20}},     //   1.0 Ap[abic] W[daef] p[fgjh] E2[cdhe] delta[ib] delta[jg] 
    {"abcd,baef,fegh,cdgh",          1.0,   4, { 6,13, 5,23}},           //   1.0 Ap[abcd] W[baef] p[fegh] E2[cdgh] 
    {"abcd,baef,fgih,cdhe,ig",       1.0,   5, { 6,13, 5,23,20}},        //   1.0 Ap[abcd] W[baef] p[fgih] E2[cdhe] delta[ig] 
    {"abic,daef,efgh,cdhg,ib",       1.0,   5, { 6,14, 5,23,20}},        //   1.0 Ap[abic] W[daef] p[efgh] E2[cdhg] delta[ib] 
    {"abic,daef,egjh,cdfh,ib,jg",    1.0,   6, { 6,14, 5,23,20,20}},     //   1.0 Ap[abic] W[daef] p[egjh] E2[cdfh] delta[ib] delta[jg] 
    {"abcd,baef,efgh,cdhg",          1.0,   4, { 6,14, 5,23}},           //   1.0 Ap[abcd] W[baef] p[efgh] E2[cdhg] 
    {"abcd,baef,egih,cdfh,ig",       1.0,   5, { 6,14, 5,23,20}},        //   1.0 Ap[abcd] W[baef] p[egih] E2[cdfh] delta[ig] 
    {"abic,dejf,aegh,cfgh,ib,jd",    2.0,   6, { 6,15, 5,23,20,21}},     //   2.0 Ap[abic] W[dejf] p[aegh] E2[cfgh] delta[ib] delta[jd] 
    {"abic,dejf,ageh,cgfh,ib,jd",   -2.0,   6, { 6,15, 5,23,20,21}},     //  -2.0 Ap[abic] W[dejf] p[ageh] E2[cgfh] delta[ib] delta[jd] 
    {"abic,dejf,aghe,cghf,ib,jd",   -2.0,   6, { 6,15, 5,23,20,21}},     //  -2.0 Ap[abic] W[dejf] p[aghe] E2[cghf] delta[ib] delta[jd] 
    {"abcd,ebif,afgh,cdgh,ie",       2.0,   5, { 6,15, 5,23,21}},        //   2.0 Ap[abcd] W[ebif] p[afgh] E2[cdgh] delta[ie] 
    {"abcd,efig,abfh,cdgh,ie",      -2.0,   5, { 6,15, 5,23,21}},        //  -2.0 Ap[abcd] W[efig] p[abfh] E2[cdgh] delta[ie] 
    {"abcd,efig,abhf,cdhg,ie",      -2.0,   5, { 6,15, 5,23,21}},        //  -2.0 Ap[abcd] W[efig] p[abhf] E2[cdhg] delta[ie] 
    {"abcd,efig,ahjf,cdgb,ie,jh",   -2.0,   6, { 6,15, 5,23,21,20}},     //  -2.0 Ap[abcd] W[efig] p[ahjf] E2[cdgb] delta[ie] delta[jh] 
    {"abic,defj,afgh,cegh,ib,jd",   -1.0,   6, { 6,16, 5,23,20,21}},     //  -1.0 Ap[abic] W[defj] p[afgh] E2[cegh] delta[ib] delta[jd] 
    {"abic,defj,ageh,cgfh,ib,jd",    1.0,   6, { 6,16, 5,23,20,21}},     //   1.0 Ap[abic] W[defj] p[ageh] E2[cgfh] delta[ib] delta[jd] 
    {"abic,defj,aghe,cghf,ib,jd",    1.0,   6, { 6,16, 5,23,20,21}},     //   1.0 Ap[abic] W[defj] p[aghe] E2[cghf] delta[ib] delta[jd] 
    {"abcd,ebfi,afgh,cdgh,ie",      -1.0,   5, { 6,16, 5,23,21}},        //  -1.0 Ap[abcd] W[ebfi] p[afgh] E2[cdgh] delta[ie] 
    {"abcd,efgi,abfh,cdgh,ie",       1.0,   5, { 6,16, 5,23,21}},        //   1.0 Ap[abcd] W[efgi] p[abfh] E2[cdgh] delta[ie] 
    {"abcd,efgi,abhf,cdhg,ie",       1.0,   5, { 6,16, 5,23,21}},        //   1.0 Ap[abcd] W[efgi] p[abhf] E2[cdhg] delta[ie] 
    {"abcd,efgi,ahjf,cdgb,ie,jh",    1.0,   6, { 6,16, 5,23,21,20}},     //   1.0 Ap[abcd] W[efgi] p[ahjf] E2[cdgb] delta[ie] delta[jh] 
    {"abic,daje,efgh,cfgh,ib,jd",    2.0,   6, { 6,17, 5,23,20,21}},     //   2.0 Ap[abic] W[daje] p[efgh] E2[cfgh] delta[ib] delta[jd] 
    {"abcd,eaif,fbgh,cdgh,ie",       2.0,   5, { 6,17, 5,23,21}},        //   2.0 Ap[abcd] W[eaif] p[fbgh] E2[cdgh] delta[ie] 
    {"abcd,eaif,fgjh,cdhb,ie,jg",    2.0,   6, { 6,17, 5,23,21,20}},     //   2.0 Ap[abcd] W[eaif] p[fgjh] E2[cdhb] delta[ie] delta[jg] 
    {"abic,daej,efgh,cfgh,ib,jd",   -1.0,   6, { 6,18, 5,23,20,21}},     //  -1.0 Ap[abic] W[daej] p[efgh] E2[cfgh] delta[ib] delta[jd] 
    {"abcd,eafi,fbgh,cdgh,ie",      -1.0,   5, { 6,18, 5,23,21}},        //  -1.0 Ap[abcd] W[eafi] p[fbgh] E2[cdgh] delta[ie] 
    {"abcd,eafi,fgjh,cdhb,ie,jg",   -1.0,   6, { 6,18, 5,23,21,20}},     //  -1.0 Ap[abcd] W[eafi] p[fgjh] E2[cdhb] delta[ie] delta[jg] 
    {"abcd,ef,aegh,cdfgbh",          1.0,   4, { 6,10, 5,25}},           //   1.0 Ap[abcd] k[ef] p[aegh] E3[cdfgbh] 
    {"abcd,ef,ageh,cdgfbh",         -1.0,   4, { 6,10, 5,25}},           //  -1.0 Ap[abcd] k[ef] p[ageh] E3[cdgfbh] 
    {"abcd,ef,aghe,cdghbf",         -1.0,   4, { 6,10, 5,25}},           //  -1.0 Ap[abcd] k[ef] p[aghe] E3[cdghbf] 
    {"abcd,ae,efgh,cdfgbh",          1.0,   4, { 6,11, 5,25}},           //   1.0 Ap[abcd] k[ae] p[efgh] E3[cdfgbh] 
    {"abjc,defg,adhi,cefhgi,jb",     1.0,   5, { 6,12, 5,25,20}},        //   1.0 Ap[abjc] W[defg] p[adhi] E3[cefhgi] delta[jb] 
    {"abjc,defg,ahdi,cehfgi,jb",    -1.0,   5, { 6,12, 5,25,20}},        //  -1.0 Ap[abjc] W[defg] p[ahdi] E3[cehfgi] delta[jb] 
    {"abjc,defg,ahid,cehigf,jb",    -1.0,   5, { 6,12, 5,25,20}},        //  -1.0 Ap[abjc] W[defg] p[ahid] E3[cehigf] delta[jb] 
    {"abcd,befg,aehi,cdghfi",        1.0,   4, { 6,12, 5,25}},           //   1.0 Ap[abcd] W[befg] p[aehi] E3[cdghfi] 
    {"abcd,befg,afhi,cdehig",        1.0,   4, { 6,12, 5,25}},           //   1.0 Ap[abcd] W[befg] p[afhi] E3[cdehig] 
    {"abcd,befg,ahei,cdhgfi",       -1.0,   4, { 6,12, 5,25}},           //  -1.0 Ap[abcd] W[befg] p[ahei] E3[cdhgfi] 
    {"abcd,befg,ahie,cdhifg",       -1.0,   4, { 6,12, 5,25}},           //  -1.0 Ap[abcd] W[befg] p[ahie] E3[cdhifg] 
    {"abcd,efgh,abei,cdfgih",       -1.0,   4, { 6,12, 5,25}},           //  -1.0 Ap[abcd] W[efgh] p[abei] E3[cdfgih] 
    {"abcd,efgh,abie,cdfigh",       -1.0,   4, { 6,12, 5,25}},           //  -1.0 Ap[abcd] W[efgh] p[abie] E3[cdfigh] 
    {"abcd,efgh,aief,cdigbh",       -1.0,   4, { 6,12, 5,25}},           //  -1.0 Ap[abcd] W[efgh] p[aief] E3[cdigbh] 
    {"abcd,efgh,aije,cdfgbh,ji",    -1.0,   5, { 6,12, 5,25,20}},        //  -1.0 Ap[abcd] W[efgh] p[aije] E3[cdfgbh] delta[ji] 
    {"abjc,daef,fghi,cdghei,jb",     1.0,   5, { 6,13, 5,25,20}},        //   1.0 Ap[abjc] W[daef] p[fghi] E3[cdghei] delta[jb] 
    {"abcd,baef,fghi,cdghei",        1.0,   4, { 6,13, 5,25}},           //   1.0 Ap[abcd] W[baef] p[fghi] E3[cdghei] 
    {"abcd,eafg,gbhi,cdehif",        1.0,   4, { 6,13, 5,25}},           //   1.0 Ap[abcd] W[eafg] p[gbhi] E3[cdehif] 
    {"abcd,eafg,gehi,cdfhbi",        1.0,   4, { 6,13, 5,25}},           //   1.0 Ap[abcd] W[eafg] p[gehi] E3[cdfhbi] 
    {"abcd,eafg,ghji,cdeibf,jh",     1.0,   5, { 6,13, 5,25,20}},        //   1.0 Ap[abcd] W[eafg] p[ghji] E3[cdeibf] delta[jh] 
    {"abjc,daef,eghi,cdgfhi,jb",     1.0,   5, { 6,14, 5,25,20}},        //   1.0 Ap[abjc] W[daef] p[eghi] E3[cdgfhi] delta[jb] 
    {"abcd,baef,eghi,cdgfhi",        1.0,   4, { 6,14, 5,25}},           //   1.0 Ap[abcd] W[baef] p[eghi] E3[cdgfhi] 
    {"abcd,eafg,fbhi,cdegih",        1.0,   4, { 6,14, 5,25}},           //   1.0 Ap[abcd] W[eafg] p[fbhi] E3[cdegih] 
    {"abcd,eafg,fghi,cdeibh",        1.0,   4, { 6,14, 5,25}},           //   1.0 Ap[abcd] W[eafg] p[fghi] E3[cdeibh] 
    {"abcd,eafg,fhji,cdegbi,jh",     1.0,   5, { 6,14, 5,25,20}},        //   1.0 Ap[abcd] W[eafg] p[fhji] E3[cdegbi] delta[jh] 
    {"abcd,efjg,afhi,cdghbi,je",     2.0,   5, { 6,15, 5,25,21}},        //   2.0 Ap[abcd] W[efjg] p[afhi] E3[cdghbi] delta[je] 
    {"abcd,efjg,ahfi,cdhgbi,je",    -2.0,   5, { 6,15, 5,25,21}},        //  -2.0 Ap[abcd] W[efjg] p[ahfi] E3[cdhgbi] delta[je] 
    {"abcd,efjg,ahif,cdhibg,je",    -2.0,   5, { 6,15, 5,25,21}},        //  -2.0 Ap[abcd] W[efjg] p[ahif] E3[cdhibg] delta[je] 
    {"abcd,efgj,aghi,cdfhbi,je",    -1.0,   5, { 6,16, 5,25,21}},        //  -1.0 Ap[abcd] W[efgj] p[aghi] E3[cdfhbi] delta[je] 
    {"abcd,efgj,ahfi,cdhgbi,je",     1.0,   5, { 6,16, 5,25,21}},        //   1.0 Ap[abcd] W[efgj] p[ahfi] E3[cdhgbi] delta[je] 
    {"abcd,efgj,ahif,cdhibg,je",     1.0,   5, { 6,16, 5,25,21}},        //   1.0 Ap[abcd] W[efgj] p[ahif] E3[cdhibg] delta[je] 
    {"abcd,eajf,fghi,cdghbi,je",     2.0,   5, { 6,17, 5,25,21}},        //   2.0 Ap[abcd] W[eajf] p[fghi] E3[cdghbi] delta[je] 
    {"abcd,eafj,fghi,cdghbi,je",    -1.0,   5, { 6,18, 5,25,21}},        //  -1.0 Ap[abcd] W[eafj] p[fghi] E3[cdghbi] delta[je] 
    {"abcd,efgh,aeij,cdfgibhj",      1.0,   4, { 6,12, 5,26}},           //   1.0 Ap[abcd] W[efgh] p[aeij] E4[cdfgibhj] 
    {"abcd,efgh,aiej,cdfigbhj",     -1.0,   4, { 6,12, 5,26}},           //  -1.0 Ap[abcd] W[efgh] p[aiej] E4[cdfigbhj] 
    {"abcd,efgh,aije,cdfijbhg",     -1.0,   4, { 6,12, 5,26}},           //  -1.0 Ap[abcd] W[efgh] p[aije] E4[cdfijbhg] 
    {"abcd,eafg,ghij,cdehibfj",      1.0,   4, { 6,13, 5,26}},           //   1.0 Ap[abcd] W[eafg] p[ghij] E4[cdehibfj] 
    {"abcd,eafg,fhij,cdehgbij",      1.0,   4, { 6,14, 5,26}},           //   1.0 Ap[abcd] W[eafg] p[fhij] E4[cdehgbij] 
	};

  FEqInfo bVec[5] = {
    {"AUVW,APQR,WR,QP,VU",   1.0,   5, { 3,19,22,20,20}}, // 1.0 b[auvw] V[apqr] E1[wr] delta[qp] delta[vu]
    {"AUVW,APQR,VWRU,QP" ,   1.0,   4, { 3,19,23,20}   }, // 1.0 b[auvw] V[apqr] E2[vwru] delta[qp]
    {"AUVW,APQR,VWQR,PU" ,   1.0,   4, { 3,19,23,20}   }, // 1.0 b[auvw] V[apqr] E2[vwqr] delta[pu]
    {"AUVW,APQR,PWRQ,VU" ,   1.0,   4, { 3,19,23,20}   }, // 1.0 b[auvw] V[apqr] E2[pwrq] delta[vu]
    {"AUVW,APQR,PVWRQU"  ,   1.0,   3, { 3,19,25}      }, // 1.0 b[auvw] V[apqr] E3[pvwrqu]
  };

  static void GetMethodInfo(FMethodInfo &Out) {
    Out = FMethodInfo();
    Out.pName = "MRLCC_AAAV";
    Out.perturberClass = "AAAV";
    Out.pSpinClass = "restricted";
    Out.pTensorDecls = &TensorDecls[0];
    Out.nTensorDecls = 27;
    Out.EqsRes = FEqSet(&EqsRes[0], 86, "MRLCC_AAAV/Res");
    Out.bVec = FEqSet(&bVec[0], 5, "MRLCC_AAAV/bVec");
  };
};
