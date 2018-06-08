namespace MRLCC_ACVV {

  FTensorDecl TensorDecls[27] = {
    /*  0*/{"R"      , "eeca"     , "", USAGE_Residual    },
    /*  1*/{"t"      , "eeca"     , "", USAGE_Amplitude   },
    /*  2*/{"T"      , "eeca"     , "", USAGE_Amplitude   },
    /*  3*/{"b"      , "eeca"     , "", USAGE_Amplitude   },
    /*  4*/{"B"      , "eeca"     , "", USAGE_Amplitude   },
    /*  5*/{"p"      , "eeca"     , "", USAGE_Amplitude   },
    /*  6*/{"Ap"     , "eeca"     , "", USAGE_Amplitude   },
    /*  7*/{"P"      , "eeca"     , "", USAGE_Amplitude   },
    /*  8*/{"AP"     , "eeca"     , "", USAGE_Amplitude   },
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
    /* 20*/{"W"      , "e"        , "", USAGE_Intermediate},
    /* 21*/{"W"      , "eeca"     , "", USAGE_Hamiltonian },
    /* 22*/{"delta"  , "cc"       , "", USAGE_Density     },
    /* 23*/{"E1"     , "aa"       , "", USAGE_Density     },
    /* 24*/{"E2"     , "aaaa"     , "", USAGE_Density     },
    /* 25*/{"S1"     , "AA"       , "", USAGE_Density     },
    /* 26*/{"S2"     , "aa"       , "", USAGE_Density     },
  };

	FEqInfo EqsRes[46] = {
    {"abcd,ef,abce,df",       -2.0,   4, { 6, 9, 5,23}},     //  -2.0 Ap[abcd] k[ef] p[abce] E1[df] 
    {"abcd,ef,bace,df",        1.0,   4, { 6, 9, 5,23}},     //   1.0 Ap[abcd] k[ef] p[bace] E1[df] 
    {"abcd,ce,abef,df",       -2.0,   4, { 6,10, 5,23}},     //  -2.0 Ap[abcd] k[ce] p[abef] E1[df] 
    {"abcd,ce,baef,df",        1.0,   4, { 6,10, 5,23}},     //   1.0 Ap[abcd] k[ce] p[baef] E1[df] 
    {"abcd,ae,becf,df",       -1.0,   4, { 6,11, 5,23}},     //  -1.0 Ap[abcd] k[ae] p[becf] E1[df] 
    {"abcd,ae,ebcf,df",        2.0,   4, { 6,11, 5,23}},     //   2.0 Ap[abcd] k[ae] p[ebcf] E1[df] 
    {"abcd,be,aecf,df",        2.0,   4, { 6,11, 5,23}},     //   2.0 Ap[abcd] k[be] p[aecf] E1[df] 
    {"abcd,be,eacf,df",       -1.0,   4, { 6,11, 5,23}},     //  -1.0 Ap[abcd] k[be] p[eacf] E1[df] 
    {"abcd,cehf,abfg,dg,he",   2.0,   5, { 6,15, 5,23,22}},  //   2.0 Ap[abcd] W[cehf] p[abfg] E1[dg] delta[he] 
    {"abcd,cehf,bafg,dg,he",  -1.0,   5, { 6,15, 5,23,22}},  //  -1.0 Ap[abcd] W[cehf] p[bafg] E1[dg] delta[he] 
    {"abcd,cefh,abfg,dg,he",  -4.0,   5, { 6,15, 5,23,22}},  //  -4.0 Ap[abcd] W[cefh] p[abfg] E1[dg] delta[he] 
    {"abcd,cefh,bafg,dg,he",   2.0,   5, { 6,15, 5,23,22}},  //   2.0 Ap[abcd] W[cefh] p[bafg] E1[dg] delta[he] 
    {"abcd,efhg,abcf,dg,he",  -4.0,   5, { 6,16, 5,23,22}},  //  -4.0 Ap[abcd] W[efhg] p[abcf] E1[dg] delta[he] 
    {"abcd,efhg,bacf,dg,he",   2.0,   5, { 6,16, 5,23,22}},  //   2.0 Ap[abcd] W[efhg] p[bacf] E1[dg] delta[he] 
    {"abcd,efgh,abcf,dg,he",   2.0,   5, { 6,17, 5,23,22}},  //   2.0 Ap[abcd] W[efgh] p[abcf] E1[dg] delta[he] 
    {"abcd,efgh,bacf,dg,he",  -1.0,   5, { 6,17, 5,23,22}},  //  -1.0 Ap[abcd] W[efgh] p[bacf] E1[dg] delta[he] 
    {"abcd,caef,bfeg,dg",      1.0,   4, { 6,18, 5,23}},     //   1.0 Ap[abcd] W[caef] p[bfeg] E1[dg] 
    {"abcd,caef,fbeg,dg",     -2.0,   4, { 6,18, 5,23}},     //  -2.0 Ap[abcd] W[caef] p[fbeg] E1[dg] 
    {"abcd,cbef,afeg,dg",     -2.0,   4, { 6,18, 5,23}},     //  -2.0 Ap[abcd] W[cbef] p[afeg] E1[dg] 
    {"abcd,cbef,faeg,dg",      1.0,   4, { 6,18, 5,23}},     //   1.0 Ap[abcd] W[cbef] p[faeg] E1[dg] 
    {"abcd,eahf,bfcg,dg,he",  -2.0,   5, { 6,18, 5,23,22}},  //  -2.0 Ap[abcd] W[eahf] p[bfcg] E1[dg] delta[he] 
    {"abcd,eahf,fbcg,dg,he",   4.0,   5, { 6,18, 5,23,22}},  //   4.0 Ap[abcd] W[eahf] p[fbcg] E1[dg] delta[he] 
    {"abcd,ebhf,afcg,dg,he",   4.0,   5, { 6,18, 5,23,22}},  //   4.0 Ap[abcd] W[ebhf] p[afcg] E1[dg] delta[he] 
    {"abcd,ebhf,facg,dg,he",  -2.0,   5, { 6,18, 5,23,22}},  //  -2.0 Ap[abcd] W[ebhf] p[facg] E1[dg] delta[he] 
    {"abcd,eafc,bfeg,dg",     -2.0,   4, { 6,19, 5,23}},     //  -2.0 Ap[abcd] W[eafc] p[bfeg] E1[dg] 
    {"abcd,eafc,fbeg,dg",      4.0,   4, { 6,19, 5,23}},     //   4.0 Ap[abcd] W[eafc] p[fbeg] E1[dg] 
    {"abcd,eafh,bfcg,dg,he",   1.0,   5, { 6,19, 5,23,22}},  //   1.0 Ap[abcd] W[eafh] p[bfcg] E1[dg] delta[he] 
    {"abcd,eafh,fbcg,dg,he",  -2.0,   5, { 6,19, 5,23,22}},  //  -2.0 Ap[abcd] W[eafh] p[fbcg] E1[dg] delta[he] 
    {"abcd,ebfc,afeg,dg",      1.0,   4, { 6,19, 5,23}},     //   1.0 Ap[abcd] W[ebfc] p[afeg] E1[dg] 
    {"abcd,ebfc,faeg,dg",     -2.0,   4, { 6,19, 5,23}},     //  -2.0 Ap[abcd] W[ebfc] p[faeg] E1[dg] 
    {"abcd,ebfh,afcg,dg,he",  -2.0,   5, { 6,19, 5,23,22}},  //  -2.0 Ap[abcd] W[ebfh] p[afcg] E1[dg] delta[he] 
    {"abcd,ebfh,facg,dg,he",   1.0,   5, { 6,19, 5,23,22}},  //   1.0 Ap[abcd] W[ebfh] p[facg] E1[dg] delta[he] 
//  {"abcd,abef,efcg,dg",      2.0,   4, { 6,20, 5,23}},     //   2.0 Ap[abcd] W[abef] p[efcg] E1[dg] 
//  {"abcd,abef,fecg,dg",     -1.0,   4, { 6,20, 5,23}},     //  -1.0 Ap[abcd] W[abef] p[fecg] E1[dg] 
    {"abcd,efgh,abce,dfgh",   -2.0,   4, { 6,12, 5,24}},     //  -2.0 Ap[abcd] W[efgh] p[abce] E2[dfgh] 
    {"abcd,efgh,bace,dfgh",    1.0,   4, { 6,12, 5,24}},     //   1.0 Ap[abcd] W[efgh] p[bace] E2[dfgh] 
    {"abcd,eafg,bgch,dehf",   -1.0,   4, { 6,13, 5,24}},     //  -1.0 Ap[abcd] W[eafg] p[bgch] E2[dehf] 
    {"abcd,eafg,gbch,dehf",    2.0,   4, { 6,13, 5,24}},     //   2.0 Ap[abcd] W[eafg] p[gbch] E2[dehf] 
    {"abcd,ebfg,agch,dehf",    2.0,   4, { 6,13, 5,24}},     //   2.0 Ap[abcd] W[ebfg] p[agch] E2[dehf] 
    {"abcd,ebfg,gach,dehf",   -1.0,   4, { 6,13, 5,24}},     //  -1.0 Ap[abcd] W[ebfg] p[gach] E2[dehf] 
    {"abcd,eafg,bfch,degh",   -1.0,   4, { 6,14, 5,24}},     //  -1.0 Ap[abcd] W[eafg] p[bfch] E2[degh] 
    {"abcd,eafg,fbch,dehg",   -1.0,   4, { 6,14, 5,24}},     //  -1.0 Ap[abcd] W[eafg] p[fbch] E2[dehg] 
    {"abcd,ebfg,afch,degh",    2.0,   4, { 6,14, 5,24}},     //   2.0 Ap[abcd] W[ebfg] p[afch] E2[degh] 
    {"abcd,ebfg,fach,degh",   -1.0,   4, { 6,14, 5,24}},     //  -1.0 Ap[abcd] W[ebfg] p[fach] E2[degh] 
    {"abcd,cefg,abfh,dehg",   -2.0,   4, { 6,16, 5,24}},     //  -2.0 Ap[abcd] W[cefg] p[abfh] E2[dehg] 
    {"abcd,cefg,bafh,dehg",    1.0,   4, { 6,16, 5,24}},     //   1.0 Ap[abcd] W[cefg] p[bafh] E2[dehg] 
    {"abcd,efgc,abeh,dfhg",    1.0,   4, { 6,17, 5,24}},     //   1.0 Ap[abcd] W[efgc] p[abeh] E2[dfhg] 
    {"abcd,efgc,baeh,dfgh",    1.0,   4, { 6,17, 5,24}},     //   1.0 Ap[abcd] W[efgc] p[baeh] E2[dfgh] 
	};

  FEqInfo bVec[2] = {
    {"ABIP,RP,ABIR",   2.0,   3, { 3,23,21}}, //  2.0 b[ABIP] E1[RP] V[ABIR]
    {"ABIP,RP,BAIR",  -1.0,   3, { 3,23,21}}, // -1.0 b[ABIP] E1[RP] V[BAIR]
  };

  int f(int i) {
    return 2*i;
  }
  FDomainDecl DomainDecls[1] = {
    {"A", "a", f}
  };

  static void GetMethodInfo(FMethodInfo &Out) {
    Out = FMethodInfo();
    Out.pName = "MRLCC_ACVV";
    Out.perturberClass = "ACVV";
    Out.pSpinClass = "restricted";
    Out.pTensorDecls = &TensorDecls[0];
    Out.nTensorDecls = 27;
    Out.pDomainDecls = &DomainDecls[0];
    Out.nDomainDecls = 1;
    Out.EqsRes = FEqSet(&EqsRes[0], 46, "MRLCC_ACVV/Res");
    Out.bVec = FEqSet(&bVec[0], 2, "MRLCC_ACVV/bVec");
  };
};
