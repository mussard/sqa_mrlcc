namespace MRLCC_AAVV {

  FTensorDecl TensorDecls[25] = {
    /*  0*/{"R"      , "eeaa"     , "", USAGE_Residual    },
    /*  1*/{"t"      , "eeaa"     , "", USAGE_Amplitude   },
    /*  2*/{"T"      , "eeaa"     , "", USAGE_Amplitude   },
    /*  3*/{"b"      , "eeaa"     , "", USAGE_Amplitude   },
    /*  4*/{"B"      , "eeaa"     , "", USAGE_Amplitude   },
    /*  5*/{"p"      , "eeaa"     , "", USAGE_Amplitude   },
    /*  6*/{"Ap"     , "eeaa"     , "", USAGE_Amplitude   },
    /*  7*/{"P"      , "eeaa"     , "", USAGE_Amplitude   },
    /*  8*/{"AP"     , "eeaa"     , "", USAGE_Amplitude   },
    /*  9*/{"k"      , "aa"       , "", USAGE_Hamiltonian },
    /* 10*/{"k"      , "ee"       , "", USAGE_Hamiltonian },
    /* 11*/{"W"      , "aaaa"     , "", USAGE_Hamiltonian },
    /* 12*/{"W"      , "aeae"     , "", USAGE_Hamiltonian },
    /* 13*/{"W"      , "aeea"     , "", USAGE_Hamiltonian },
    /* 14*/{"W"      , "caca"     , "", USAGE_Hamiltonian },
    /* 15*/{"W"      , "caac"     , "", USAGE_Hamiltonian },
    /* 16*/{"W"      , "cece"     , "", USAGE_Hamiltonian },
    /* 17*/{"W"      , "ceec"     , "", USAGE_Hamiltonian },
    /* 18*/{"W"      , "e"        , "", USAGE_Intermediate},
    /* 19*/{"W"      , "eeaa"     , "", USAGE_Hamiltonian },
    /* 20*/{"delta"  , "cc"       , "", USAGE_Density     },
    /* 21*/{"E2"     , "aaaa"     , "", USAGE_Density     },
    /* 22*/{"S1"     , "aaaa"     , "", USAGE_Density     },
    /* 23*/{"S2"     , "aaaa"     , "", USAGE_Density     },
    /* 24*/{"E3"     , "aaaaaa"   , "", USAGE_Density     },
  };

	FEqInfo EqsRes[10] = {
    {"abcd,ef,abeg,cdfg",       -4.0,   4, { 6, 9, 5,21}},     //  -4.0 Ap[abcd] k[ef] p[abeg] E2[cdfg] 
    {"abcd,ae,befg,cdgf",        4.0,   4, { 6,10, 5,21}},     //   4.0 Ap[abcd] k[ae] p[befg] E2[cdgf] 
    {"abcd,efgh,abef,cdgh",     -2.0,   4, { 6,11, 5,21}},     //  -2.0 Ap[abcd] W[efgh] p[abef] E2[cdgh] 
    {"abcd,efig,abfh,cdgh,ie",  -8.0,   5, { 6,14, 5,21,20}},  //  -8.0 Ap[abcd] W[efig] p[abfh] E2[cdgh] delta[ie] 
    {"abcd,efgi,abfh,cdgh,ie",   4.0,   5, { 6,15, 5,21,20}},  //   4.0 Ap[abcd] W[efgi] p[abfh] E2[cdgh] delta[ie] 
    {"abcd,eaif,bfgh,cdhg,ie",   8.0,   5, { 6,16, 5,21,20}},  //   8.0 Ap[abcd] W[eaif] p[bfgh] E2[cdhg] delta[ie] 
    {"abcd,eafi,bfgh,cdhg,ie",  -4.0,   5, { 6,17, 5,21,20}},  //  -4.0 Ap[abcd] W[eafi] p[bfgh] E2[cdhg] delta[ie] 
//  {"abcd,abef,efgh,cdgh",      2.0,   4, { 6,18, 5,21}},     //   2.0 Ap[abcd] W[abef] p[efgh] E2[cdgh] 
    {"abcd,efgh,abei,cdfgih",   -4.0,   4, { 6,11, 5,24}},     //  -4.0 Ap[abcd] W[efgh] p[abei] E3[cdfgih] 
    {"abcd,eafg,bghi,cdeihf",    4.0,   4, { 6,12, 5,24}},     //   4.0 Ap[abcd] W[eafg] p[bghi] E3[cdeihf] 
    {"abcd,eafg,bfhi,cdeghi",    4.0,   4, { 6,13, 5,24}},     //   4.0 Ap[abcd] W[eafg] p[bfhi] E3[cdeghi] 
	};

  FEqInfo bVec[4] = {
    {"ABRS,ABPQ,RSPQ",  0.25,   3, { 3,19,21}}, // 0.5 b[ABRS] V[ABPQ] E2[RSPQ]
    {"BARS,ABPQ,RSQP",  0.25,   3, { 3,19,21}}, // 0.5 b[BARS] V[ABPQ] E2[RSQP]
    {"ABRS,BAPQ,RSQP",  0.25,   3, { 3,19,21}}, // 0.5 b[ABRS] V[BAPQ] E2[RSQP]
    {"BARS,BAPQ,RSPQ",  0.25,   3, { 3,19,21}}, // 0.5 b[BARS] V[BAPQ] E2[RSPQ]
  };

  static void GetMethodInfo(FMethodInfo &Out) {
    Out = FMethodInfo();
    Out.pName = "MRLCC_AAVV";
    Out.perturberClass = "AAVV";
    Out.pSpinClass = "restricted";
    Out.pTensorDecls = &TensorDecls[0];
    Out.nTensorDecls = 25;
    Out.EqsRes = FEqSet(&EqsRes[0], 10, "MRLCC_AAVV/Res");
    Out.bVec = FEqSet(&bVec[0], 4, "MRLCC_AAVV/bVec");
  };
};
