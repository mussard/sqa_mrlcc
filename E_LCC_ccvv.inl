namespace MRLCC_CCVV {

  FTensorDecl TensorDecls[22] = {
    /*  0*/{"R"      , "eecc"     , "", USAGE_Residual    },
    /*  1*/{"t"      , "eecc"     , "", USAGE_Amplitude   },
    /*  2*/{"T"      , "eecc"     , "", USAGE_Amplitude   },
    /*  3*/{"b"      , "eecc"     , "", USAGE_Amplitude   },
    /*  4*/{"B"      , "eecc"     , "", USAGE_Amplitude   },
    /*  5*/{"p"      , "eecc"     , "", USAGE_Amplitude   },
    /*  6*/{"Ap"     , "eecc"     , "", USAGE_Amplitude   },
    /*  7*/{"P"      , "eecc"     , "", USAGE_Amplitude   },
    /*  8*/{"AP"     , "eecc"     , "", USAGE_Amplitude   },
    /*  9*/{"k"      , "cc"       , "", USAGE_Hamiltonian },
    /* 10*/{"k"      , "ee"       , "", USAGE_Hamiltonian },
    /* 11*/{"W"      , "aeae"     , "", USAGE_Hamiltonian },
    /* 12*/{"W"      , "aeea"     , "", USAGE_Hamiltonian },
    /* 13*/{"W"      , "cccc"     , "", USAGE_Hamiltonian },
    /* 14*/{"W"      , "caca"     , "", USAGE_Hamiltonian },
    /* 15*/{"W"      , "caac"     , "", USAGE_Hamiltonian },
    /* 16*/{"W"      , "cece"     , "", USAGE_Hamiltonian },
    /* 17*/{"W"      , "ceec"     , "", USAGE_Hamiltonian },
    /* 18*/{"W"      , "e"        , "", USAGE_Intermediate},
    /* 19*/{"W"      , "eecc"     , "", USAGE_Hamiltonian },
    /* 20*/{"delta"  , "cc"       , "", USAGE_Density     },
    /* 21*/{"E1"     , "aa"       , "", USAGE_Density     },
  };

	FEqInfo EqsRes[30] = {
    {"abcd,ce,abde",        8.0,   3, { 6, 9, 5}},     //   8.0 Ap[abcd] k[ce] p[abde] 
    {"abcd,ce,abed",      -16.0,   3, { 6, 9, 5}},     // -16.0 Ap[abcd] k[ce] p[abed] 
    {"abcd,ae,becd",       -8.0,   3, { 6,10, 5}},     //  -8.0 Ap[abcd] k[ae] p[becd] 
    {"abcd,ae,bedc",       16.0,   3, { 6,10, 5}},     //  16.0 Ap[abcd] k[ae] p[bedc] 
    {"abcd,cdef,abef",      8.0,   3, { 6,13, 5}},     //   8.0 Ap[abcd] W[cdef] p[abef] 
    {"abcd,cdef,abfe",     -4.0,   3, { 6,13, 5}},     //  -4.0 Ap[abcd] W[cdef] p[abfe] 
    {"abcd,cegf,abdf,ge",  -8.0,   4, { 6,13, 5,20}},  //  -8.0 Ap[abcd] W[cegf] p[abdf] delta[ge] 
    {"abcd,cegf,abfd,ge",  16.0,   4, { 6,13, 5,20}},  //  16.0 Ap[abcd] W[cegf] p[abfd] delta[ge] 
    {"abcd,cefg,abdf,ge",  16.0,   4, { 6,13, 5,20}},  //  16.0 Ap[abcd] W[cefg] p[abdf] delta[ge] 
    {"abcd,cefg,abfd,ge", -32.0,   4, { 6,13, 5,20}},  // -32.0 Ap[abcd] W[cefg] p[abfd] delta[ge] 
    {"abcd,caef,bfde",    -16.0,   3, { 6,16, 5}},     // -16.0 Ap[abcd] W[caef] p[bfde] 
    {"abcd,caef,bfed",      8.0,   3, { 6,16, 5}},     //   8.0 Ap[abcd] W[caef] p[bfed] 
    {"abcd,cbef,afde",      8.0,   3, { 6,16, 5}},     //   8.0 Ap[abcd] W[cbef] p[afde] 
    {"abcd,cbef,afed",    -16.0,   3, { 6,16, 5}},     // -16.0 Ap[abcd] W[cbef] p[afed] 
    {"abcd,eagf,bfcd,ge", -16.0,   4, { 6,16, 5,20}},  // -16.0 Ap[abcd] W[eagf] p[bfcd] delta[ge] 
    {"abcd,eagf,bfdc,ge",  32.0,   4, { 6,16, 5,20}},  //  32.0 Ap[abcd] W[eagf] p[bfdc] delta[ge] 
    {"abcd,eafc,bfde",     32.0,   3, { 6,17, 5}},     //  32.0 Ap[abcd] W[eafc] p[bfde] 
    {"abcd,eafc,bfed",    -16.0,   3, { 6,17, 5}},     // -16.0 Ap[abcd] W[eafc] p[bfed] 
    {"abcd,eafd,bfce",    -16.0,   3, { 6,17, 5}},     // -16.0 Ap[abcd] W[eafd] p[bfce] 
    {"abcd,eafd,bfec",      8.0,   3, { 6,17, 5}},     //   8.0 Ap[abcd] W[eafd] p[bfec] 
    {"abcd,eafg,bfcd,ge",   8.0,   4, { 6,17, 5,20}},  //   8.0 Ap[abcd] W[eafg] p[bfcd] delta[ge] 
    {"abcd,eafg,bfdc,ge", -16.0,   4, { 6,17, 5,20}},  // -16.0 Ap[abcd] W[eafg] p[bfdc] delta[ge] 
//  {"abcd,abef,efcd",      8.0,   3, { 6,18, 5}},     //   8.0 Ap[abcd] W[abef] p[efcd] 
//  {"abcd,abef,efdc",     -4.0,   3, { 6,18, 5}},     //  -4.0 Ap[abcd] W[abef] p[efdc] 
    {"abcd,eafg,bgcd,ef",  -8.0,   4, { 6,11, 5,21}},  //  -8.0 Ap[abcd] W[eafg] p[bgcd] E1[ef] 
    {"abcd,eafg,bgdc,ef",  16.0,   4, { 6,11, 5,21}},  //  16.0 Ap[abcd] W[eafg] p[bgdc] E1[ef] 
    {"abcd,eafg,bfcd,eg",   4.0,   4, { 6,12, 5,21}},  //   4.0 Ap[abcd] W[eafg] p[bfcd] E1[eg] 
    {"abcd,eafg,bfdc,eg",  -8.0,   4, { 6,12, 5,21}},  //  -8.0 Ap[abcd] W[eafg] p[bfdc] E1[eg] 
    {"abcd,cefg,abdf,eg",   8.0,   4, { 6,14, 5,21}},  //   8.0 Ap[abcd] W[cefg] p[abdf] E1[eg] 
    {"abcd,cefg,abfd,eg", -16.0,   4, { 6,14, 5,21}},  // -16.0 Ap[abcd] W[cefg] p[abfd] E1[eg] 
    {"abcd,efgc,abde,fg",  -4.0,   4, { 6,15, 5,21}},  //  -4.0 Ap[abcd] W[efgc] p[abde] E1[fg] 
    {"abcd,efgc,abed,fg",   8.0,   4, { 6,15, 5,21}},  //   8.0 Ap[abcd] W[efgc] p[abed] E1[fg] 
	};

  FEqInfo bVec[4] = {
    {"CDKL,LM,CDKM",   2.0,   3, { 3,20,19}}, //  2.0 b[CDKL] deltac[LM] V[CDKM]
    {"CDKL,LM,DCKM",  -1.0,   3, { 3,20,19}}, // -1.0 b[CDKL] deltac[LM] V[DCKM]
    {"CDKL,LM,CDMK",  -1.0,   3, { 3,20,19}}, // -1.0 b[CDKL] deltac[LM] V[CDMK]
    {"CDKL,LM,DCMK",   2.0,   3, { 3,20,19}}, //  2.0 b[CDKL] deltac[LM] V[DCMK]
  };

  static void GetMethodInfo(FMethodInfo &Out) {
    Out = FMethodInfo();
    Out.pName = "MRLCC_CCVV";
    Out.perturberClass = "CCVV";
    Out.pSpinClass = "restricted";
    Out.pTensorDecls = &TensorDecls[0];
    Out.nTensorDecls = 22;
    Out.EqsRes = FEqSet(&EqsRes[0], 30, "MRLCC_CCVV/Res");
    Out.bVec = FEqSet(&bVec[0], 4, "MRLCC_CCVV/bVec");
  };
};
