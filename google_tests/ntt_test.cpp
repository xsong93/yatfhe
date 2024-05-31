//
// Created by Xintong Song on 2024/1/30.
//
#include "gtest/gtest.h"
#include "yatfhe/ntt.h"
#include "yatfhe/ntt64.h"
#include "yatfhe/polynomial.h"
#include "yatfhe/numeric_functions.h"
#include "yautil/time_counter.h"
#include "yautil/tool.h"

std::vector<int> extractValues(const std::string& input) {
    std::vector<int> result;
    std::istringstream iss(input);
    std::string token;

    while (std::getline(iss, token, ' ')) {
        size_t pos = token.find(':');
        if (pos != std::string::npos) {
            std::string valueStr = token.substr(pos + 1);
            int value = std::stoi(valueStr);
            result.push_back(value);
        }
    }

    return result;
}

TEST(NttTest, NttIntt64Test) {
    string a = "[0:-2 1:7 2:4 3:8 4:-7 5:4 6:6 7:0 8:-1 9:-4 10:-4 11:-1 12:-16 13:7 14:-1 15:-2 16:0 17:-2 18:-4 19:-7 20:-1 21:1 22:3 23:-8 24:3 25:2 26:7 27:-4 28:6 29:-7 30:-4 31:1 32:7 33:5 34:6 35:-7 36:3 37:6 38:7 39:-16 40:-3 41:3 42:3 43:7 44:8 45:-2 46:-4 47:-2 48:-3 49:-1 50:-1 51:-3 52:5 53:-4 54:2 55:-7 56:1 57:-8 58:5 59:4 60:-16 61:-6 62:-3 63:5 64:4 65:-4 66:-1 67:-5 68:0 69:8 70:4 71:-7 72:5 73:8 74:7 75:3 76:-2 77:-5 78:3 79:16 80:3 81:7 82:6 83:6 84:-5 85:-2 86:8 87:-3 88:-6 89:-7 90:-1 91:-6 92:5 93:8 94:4 95:4 96:-1 97:1 98:0 99:4 100:-7 101:-7 102:5 103:-7 104:-5 105:-1 106:-3 107:-3 108:-1 109:-3 110:-8 111:0 112:-4 113:-2 114:6 115:2 116:-4 117:-16 118:4 119:-1 120:2 121:-7 122:5 123:-2 124:-8 125:3 126:3 127:1 128:8 129:7 130:-4 131:2 132:-5 133:-6 134:-3 135:-2 136:3 137:2 138:-1 139:3 140:7 141:-5 142:-6 143:7 144:-3 145:-3 146:4 147:6 148:-1 149:-2 150:-2 151:-5 152:3 153:-2 154:4 155:-6 156:-8 157:-5 158:7 159:-4 160:-6 161:6 162:-7 163:2 164:-6 165:7 166:-3 167:-1 168:2 169:3 170:-3 171:-5 172:8 173:4 174:7 175:2 176:-5 177:5 178:-4 179:-2 180:-7 181:-6 182:0 183:-2 184:0 185:7 186:0 187:5 188:1 189:-2 190:-6 191:1 192:-1 193:5 194:1 195:-1 196:4 197:6 198:6 199:8 200:4 201:1 202:-2 203:-1 204:-7 205:1 206:5 207:5 208:1 209:2 210:4 211:8 212:3 213:6 214:-4 215:-8 216:-4 217:1 218:4 219:-1 220:1 221:7 222:5 223:2 224:4 225:-8 226:7 227:1 228:7 229:-3 230:3 231:-6 232:-6 233:1 234:7 235:-5 236:4 237:8 238:1 239:-3 240:-6 241:6 242:2 243:3 244:6 245:-1 246:7 247:16 248:-3 249:5 250:-5 251:-4 252:6 253:5 254:6 255:4 256:-3 257:-2 258:-1 259:-4 260:1 261:-3 262:-7 263:-5 264:-5 265:2 266:3 267:0 268:-6 269:3 270:8 271:5 272:-5 273:-1 274:-16 275:5 276:-3 277:1 278:-1 279:-2 280:1 281:6 282:-4 283:-4 284:-8 285:-5 286:4 287:-8 288:6 289:4 290:7 291:7 292:-7 293:16 294:5 295:-3 296:0 297:2 298:3 299:3 300:-1 301:7 302:6 303:-2 304:-8 305:-5 306:7 307:0 308:-3 309:2 310:-5 311:4 312:-2 313:4 314:8 315:-4 316:5 317:-8 318:0 319:5 320:-4 321:-16 322:0 323:3 324:-7 325:2 326:5 327:-4 328:-6 329:6 330:-1 331:2 332:-6 333:8 334:4 335:-3 336:-6 337:0 338:-6 339:3 340:-6 341:-8 342:-4 343:5 344:-5 345:4 346:-4 347:-3 348:-4 349:-6 350:-5 351:-3 352:-16 353:-1 354:-6 355:4 356:6 357:-4 358:-4 359:-4 360:-16 361:-8 362:-2 363:7 364:-2 365:-4 366:-1 367:-8 368:7 369:2 370:-16 371:8 372:8 373:2 374:-6 375:3 376:-7 377:-2 378:2 379:-5 380:-5 381:-7 382:-6 383:1 384:-4 385:-7 386:4 387:2 388:7 389:4 390:-2 391:3 392:-3 393:-3 394:6 395:16 396:-3 397:-6 398:-1 399:-7 400:-1 401:8 402:7 403:-1 404:-4 405:-1 406:-2 407:-7 408:-5 409:-16 410:-2 411:-5 412:6 413:-3 414:-1 415:0 416:7 417:-8 418:-3 419:-8 420:4 421:7 422:6 423:-5 424:-8 425:4 426:4 427:8 428:1 429:1 430:-5 431:7 432:3 433:8 434:-6 435:5 436:-7 437:-7 438:5 439:-3 440:0 441:0 442:5 443:4 444:1 445:-8 446:8 447:-4 448:7 449:-3 450:5 451:-6 452:8 453:-8 454:-5 455:-2 456:4 457:6 458:1 459:1 460:4 461:-5 462:4 463:2 464:6 465:-4 466:4 467:4 468:-3 469:2 470:6 471:-1 472:-8 473:1 474:4 475:5 476:-2 477:4 478:-3 479:-1 480:6 481:3 482:6 483:-6 484:-6 485:-1 486:16 487:2 488:-2 489:8 490:-2 491:-16 492:-3 493:-6 494:-1 495:-2 496:2 497:-6 498:-1 499:0 500:3 501:-2 502:-6 503:5 504:2 505:3 506:5 507:-6 508:-1 509:7 510:-1 511:5 512:-3 513:-7 514:3 515:16 516:1 517:2 518:-5 519:1 520:-5 521:2 522:0 523:-3 524:7 525:-8 526:4 527:-2 528:-5 529:-5 530:7 531:-2 532:-2 533:4 534:-5 535:5 536:-2 537:4 538:-4 539:-2 540:-7 541:2 542:-6 543:-3 544:7 545:-7 546:6 547:-7 548:6 549:-7 550:-4 551:-7 552:-5 553:-2 554:-6 555:2 556:-8 557:3 558:2 559:3 560:-2 561:2 562:6 563:16 564:-8 565:-4 566:3 567:-3 568:-3 569:-2 570:-16 571:-6 572:16 573:-5 574:7 575:-2 576:0 577:-6 578:3 579:4 580:5 581:-2 582:-3 583:2 584:3 585:5 586:2 587:-1 588:2 589:5 590:-6 591:2 592:-5 593:5 594:-16 595:8 596:5 597:7 598:0 599:8 600:-3 601:3 602:-3 603:-2 604:0 605:-7 606:2 607:1 608:3 609:-2 610:5 611:2 612:-2 613:-8 614:6 615:3 616:2 617:-3 618:-2 619:6 620:-2 621:4 622:-2 623:1 624:2 625:-6 626:-16 627:3 628:8 629:0 630:5 631:-4 632:-7 633:-8 634:8 635:3 636:3 637:3 638:0 639:-6 640:6 641:8 642:-7 643:-2 644:1 645:3 646:-4 647:16 648:6 649:6 650:4 651:-6 652:-4 653:-8 654:-6 655:-16 656:-6 657:7 658:-6 659:-3 660:5 661:2 662:6 663:1 664:6 665:-6 666:-5 667:7 668:2 669:7 670:-1 671:-3 672:2 673:2 674:-5 675:4 676:-6 677:3 678:-3 679:0 680:6 681:3 682:7 683:6 684:4 685:-6 686:8 687:5 688:5 689:7 690:-5 691:5 692:3 693:-3 694:-1 695:-1 696:6 697:-1 698:7 699:-5 700:-5 701:7 702:-3 703:-2 704:1 705:-7 706:-1 707:-6 708:-2 709:-5 710:8 711:-1 712:3 713:0 714:3 715:-2 716:3 717:1 718:5 719:7 720:4 721:-6 722:1 723:1 724:-5 725:4 726:-5 727:-4 728:3 729:8 730:3 731:5 732:2 733:-6 734:6 735:1 736:7 737:-6 738:5 739:1 740:5 741:6 742:3 743:7 744:1 745:3 746:-1 747:3 748:-4 749:-5 750:-2 751:-5 752:-7 753:-7 754:-6 755:4 756:-2 757:-5 758:4 759:4 760:-16 761:-8 762:1 763:0 764:5 765:3 766:-6 767:-6 768:-8 769:7 770:6 771:3 772:0 773:-3 774:2 775:5 776:2 777:6 778:-1 779:6 780:4 781:-6 782:-8 783:-16 784:3 785:-1 786:-7 787:-6 788:4 789:-1 790:-5 791:3 792:3 793:-6 794:6 795:-6 796:-2 797:-3 798:-4 799:-6 800:8 801:5 802:8 803:3 804:-1 805:-3 806:-5 807:-5 808:0 809:4 810:7 811:-3 812:-7 813:-4 814:5 815:-1 816:7 817:-2 818:-2 819:-16 820:-4 821:-8 822:4 823:16 824:3 825:4 826:-3 827:2 828:-8 829:1 830:2 831:-8 832:7 833:-7 834:-3 835:0 836:1 837:3 838:1 839:3 840:6 841:-6 842:-1 843:1 844:16 845:8 846:5 847:7 848:-4 849:-8 850:1 851:-7 852:3 853:4 854:8 855:-16 856:3 857:3 858:-5 859:-7 860:4 861:-6 862:-3 863:2 864:-1 865:-6 866:-7 867:6 868:3 869:-6 870:-7 871:4 872:0 873:0 874:5 875:-16 876:7 877:2 878:-5 879:-6 880:1 881:8 882:1 883:3 884:5 885:-8 886:2 887:7 888:0 889:2 890:2 891:4 892:-6 893:-6 894:7 895:2 896:-4 897:4 898:7 899:-6 900:7 901:0 902:-2 903:0 904:1 905:8 906:4 907:7 908:4 909:-5 910:-3 911:7 912:-3 913:1 914:-5 915:-4 916:-2 917:-8 918:-7 919:4 920:5 921:-3 922:6 923:3 924:4 925:3 926:3 927:7 928:-3 929:-6 930:-6 931:1 932:2 933:4 934:1 935:-7 936:3 937:16 938:-2 939:-1 940:7 941:3 942:-1 943:6 944:-3 945:2 946:6 947:2 948:-3 949:8 950:-7 951:-8 952:-8 953:6 954:-8 955:2 956:8 957:-1 958:2 959:8 960:7 961:-6 962:-1 963:-5 964:-1 965:-5 966:16 967:16 968:-2 969:5 970:-4 971:4 972:0 973:5 974:2 975:-5 976:-5 977:8 978:4 979:-6 980:4 981:2 982:3 983:-4 984:4 985:7 986:-4 987:8 988:-5 989:1 990:6 991:7 992:-4 993:16 994:8 995:5 996:-6 997:4 998:7 999:7 1000:16 1001:4 1002:-6 1003:8 1004:-16 1005:-7 1006:-7 1007:3 1008:2 1009:8 1010:8 1011:5 1012:-4 1013:6 1014:-4 1015:5 1016:1 1017:-4 1018:8 1019:3 1020:6 1021:-3 1022:1 1023:-8 ]";
    std::vector<int> aPrime = extractValues(a);
    printArray(aPrime, "a");
    const int N = 1024;
    initGlobalParamsNtt64(N);
    LagrangePolynomial resNtt{N};
    TorusPolynomial a1{N};
    TorusPolynomial resIntt{N};
    initCoeffsViaUniformDistribution(a1.coeffs);
//    for (auto i = 0; i < N; i++) {
//        a1.coeffs[i] = -1;
//    }
    applyNtt(resNtt, a1);
    applyIntt(resIntt, resNtt);
    printArray(resNtt.coeffs, "resNtt");
    printArray(a1.coeffs, "a0");
    printArray(resIntt.coeffs, "a1");
    for (auto i = 0; i < a1.N; i++) {
        if (a1.coeffs[i] != resIntt.coeffs[i]) {
            printf("NE at index %d\n",i);
        }
        ASSERT_EQ(a1.coeffs[i], resIntt.coeffs[i]);
    }
    printBanner("NttIntt64Test");
}

TEST(NttTest, NttAddConstantTest) {
    const int N = 1024;
    initGlobalParamsNtt64(N);
    LagrangePolynomial a{N};
    LagrangePolynomial b{N};
    LagrangePolynomial resNtt{N};
    IntPolynomial poly{N};
    IntPolynomial c{N};
    IntPolynomial res{N};
    for (int i = 0; i < a.N; i++) {
        poly.coeffs[i] = i;
    }
    c.coeffs[0] = 77;
    applyNtt(a, poly);
    applyNtt(b, c);
    printArray(c.coeffs, "cOri");
    printArray(b.coeffs, "bNtt");
    for (int i = 0; i < a.N; i++) {
        resNtt.coeffs[i] = modAdd(a.coeffs[i], b.coeffs[i]);
    }
    applyIntt(res, resNtt);
    printArray(res.coeffs, "res");
}

TEST(NttTest, NttSamePolyTest) {
    COUNT_TIME("init timer", cout << endl;)
    const int N = 1024;
    initGlobalParamsNtt64(N);
    LagrangePolynomial a{N};
    LagrangePolynomial b{N};
    LagrangePolynomial tmpMul{N};
    LagrangePolynomial tmpAdd{N};
    LagrangePolynomial tmpSub{N};

    TorusPolynomial poly0{N};
    TorusPolynomial poly2{N};
    TorusPolynomial resMul{N};
    TorusPolynomial resAdd{N};
    TorusPolynomial resSub{N};
    TorusPolynomial navMul{N};
    TorusPolynomial navAdd{N};
    TorusPolynomial navSub{N};

    for (auto j = 0; j < N; j++) {
        poly0.coeffs[j] = genIntUniformDist(IntMin ,IntMax);
        poly2.coeffs[j] = genIntUniformDist(0 ,1);
    }
    printArray(poly0.coeffs, "poly0");
    printArray(poly2.coeffs, "poly2");

    COUNT_TIME("NTT_MULT", {
        applyNtt(a, poly0);
        applyNtt(b, poly2);
        for (int i = 0; i < a.N; i++) {
           tmpMul.coeffs[i] = modMULT64(a.coeffs[i], b.coeffs[i]);
        }
        applyIntt(resMul, tmpMul);})
    COUNT_TIME("NAIVE_MULT",
        polynomialMulNaive(navMul, poly0, poly2);)

    for (int i = 0; i < a.N; i++) {
        tmpAdd.coeffs[i] = modAdd(a.coeffs[i], b.coeffs[i]);
        tmpSub.coeffs[i] = modSub(a.coeffs[i], b.coeffs[i]);
    }
    applyIntt(resAdd, tmpAdd);
    applyIntt(resSub, tmpSub);

    polynomialAdd(navAdd, poly0, poly2);
    polynomialSub(navSub, poly0, poly2);

    for (int i = 0; i < navMul.N; i++) {
        EXPECT_EQ(resMul.coeffs[i], navMul.coeffs[i]);
        EXPECT_EQ(resAdd.coeffs[i], navAdd.coeffs[i]);
        EXPECT_EQ(resSub.coeffs[i], navSub.coeffs[i]);
    }
    printBanner("NttSamePoly");
}

TEST(NttTest, ConvolutionTest) {
    COUNT_TIME("init timer", cout << endl;)
    const int N = 512;
    const int k = 5;
    initGlobalParamsNtt64(N);

    vector<LagrangePolynomial> a(k, LagrangePolynomial(N));
    vector<LagrangePolynomial> b(k, LagrangePolynomial(N));
    LagrangePolynomial tmpMul{N};

    vector<TorusPolynomial> poly0(k, TorusPolynomial(N));
    vector<TorusPolynomial> poly2(k, TorusPolynomial(N));
    TorusPolynomial resMul{N};
    TorusPolynomial navMul{N};
    for (auto i = 0 ; i < k; i++) {
        for (auto j = 0; j < N; j++) {
            poly0[i].coeffs[j] = genIntUniformDist(IntMin ,IntMax);
            poly2[i].coeffs[j] = genIntUniformDist(0 ,1);
        }
        printArray(poly0[i].coeffs, "poly0" + to_string(i));
        printArray(poly2[i].coeffs, "poly2" + to_string(i));
    }

    COUNT_TIME("NTT_MULT", {
        for (auto i = 0 ; i < k; i++) {
            applyNtt(a[i], poly0[i]);
            applyNtt(b[i], poly2[i]);
        }
        calModularInnerProductNtt(tmpMul, a, b);
        applyIntt(resMul, tmpMul);})
    COUNT_TIME("NAIVE_MULT",
        for (auto i = 0 ; i < k; i++) {
            polynomialMulAccNaive(navMul, poly0[i], poly2[i]);
        })
    printArray(resMul.coeffs, "resMul");
    printArray(navMul.coeffs, "navMul");

    for (int i = 0; i < navMul.N; i++) {
        EXPECT_EQ(resMul.coeffs[i], navMul.coeffs[i]);
    }
    printBanner("NttSamePoly");
}


//
//TEST(NttRotTest, NttRotTest) {
//    const int N = 1024;
//    LagrangePolynomial a1Ntt{N};
//    LagrangePolynomial a2Ntt{N};
//    TorusPolynomial a1{N};
//    TorusPolynomial a2{N};
//    TorusPolynomial a1Intt{N};
//    TorusPolynomial a2Intt{N};
//    initCoeffsViaUniformDistribution(a1.coeffs);
//    for (auto i = 0; i < N; i++) {
//        a1.coeffs[i] = i + 1;
//    }
//    torusPolynomialRotate(a2, 3, a1);
//    printArray(a1.coeffs, "a1");
//    printArray(a2.coeffs, "a2");
//
//    applyNtt(a1Ntt, a1);
//    applyNtt(a2Ntt, a2);
//    printArray(a1Ntt.coeffs, "a1Ntt");
//    printArray(a2Ntt.coeffs, "a2Ntt");
//
//    LagrangePolynomial tmp{N};
//
//
//    //todo: rot
//
//    LagrangePolynomial sub{N};
//    lagrangePolynomialSub(sub, a2Ntt, a1Ntt);
//
//    applyIntt(a1Intt, sub);
//    printArray(a1Intt.coeffs, "a1Intt");
//
//    printBanner("NttRotTest");
//}

TEST(NttTest, bit_rev_test) {
    int N = 32;
    std::vector<NttType> vec(N);
    for (int i = 0; i < N; i++) {
        vec[i] = i;
    }
    for (int j = 0; j < N; j++) {
        std::cout << j << ":" << vec[j] << " ";
    }
    std::cout<<std::endl;
    bitRevShuffle(vec);
    for (int k = 0; k < N; k++) {
        std::cout << k << ":" << vec[k] << " ";
    }
    std::cout<<std::endl;

}