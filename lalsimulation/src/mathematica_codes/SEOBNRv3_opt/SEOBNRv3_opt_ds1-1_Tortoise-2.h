REAL8 tmp1=x->data[0]*x->data[0];
REAL8 tmp2=x->data[1]*x->data[1];
REAL8 tmp3=x->data[2]*x->data[2];
REAL8 tmp4=tmp1+tmp2+tmp3;
REAL8 tmp9=s1Vec->data[1]+s2Vec->data[1];
REAL8 tmp7=s1Vec->data[0]+s2Vec->data[0];
REAL8 tmp8=tmp7*tmp7;
REAL8 tmp10=tmp9*tmp9;
REAL8 tmp11=s1Vec->data[2]+s2Vec->data[2];
REAL8 tmp12=tmp11*tmp11;
REAL8 tmp13=tmp10+tmp12+tmp8;
REAL8 tmp15=1./sqrt(tmp13);
REAL8 tmp16=1./sqrt(tmp4);
REAL8 tmp5=(1.0/(tmp4*tmp4));
REAL8 tmp43=1/tmp4;
REAL8 tmp46=coeffs->KK*eta;
REAL8 tmp47=-1.+tmp46;
REAL8 tmp35=pow(tmp4,-2.5);
REAL8 tmp62=tmp13*tmp13;
REAL8 tmp41=(1.0/sqrt(tmp4*tmp4*tmp4));
REAL8 tmp48=(1.0/(tmp47*tmp47));
REAL8 tmp49=1.*tmp48;
REAL8 tmp50=1.*tmp13*tmp43;
REAL8 tmp51=1/tmp47;
REAL8 tmp52=2.*tmp16*tmp51;
REAL8 tmp53=tmp49+tmp50+tmp52;
REAL8 tmp86=(1.0/sqrt(tmp13*tmp13*tmp13));
REAL8 tmp54=tmp15*tmp16*tmp7*x->data[0];
REAL8 tmp55=tmp15*tmp16*tmp9*x->data[1];
REAL8 tmp56=tmp11*tmp15*tmp16*x->data[2];
REAL8 tmp57=tmp54+tmp55+tmp56;
REAL8 tmp76=-1.+m1PlusEtaKK;
REAL8 tmp61=c1k5*tmp13;
REAL8 tmp63=c2k5*tmp62;
REAL8 tmp64=c0k5+tmp61+tmp63;
REAL8 tmp65=1.*tmp35*tmp64;
REAL8 tmp66=c1k4*tmp13;
REAL8 tmp67=c2k4*tmp62;
REAL8 tmp68=c0k4+tmp66+tmp67;
REAL8 tmp69=1.*tmp5*tmp68;
REAL8 tmp70=c1k3*tmp13;
REAL8 tmp71=c0k3+tmp70;
REAL8 tmp72=1.*tmp41*tmp71;
REAL8 tmp73=c1k2*tmp13;
REAL8 tmp74=c0k2+tmp73;
REAL8 tmp75=1.*tmp43*tmp74;
REAL8 tmp77=coeffs->KK*tmp76;
REAL8 tmp78=coeffs->KK+tmp77;
REAL8 tmp79=-2.*m1PlusEtaKK*tmp16*tmp78;
REAL8 tmp80=1.*tmp16;
REAL8 tmp81=log(tmp80);
REAL8 tmp82=1.*coeffs->k5l*tmp35*tmp81;
REAL8 tmp83=1.+tmp65+tmp69+tmp72+tmp75+tmp79+tmp82;
REAL8 tmp58=tmp57*tmp57;
REAL8 tmp59=-tmp58;
REAL8 tmp60=1.+tmp59;
REAL8 tmp92=coeffs->KK*eta*tmp76;
REAL8 tmp93=log(tmp83);
REAL8 tmp94=eta*tmp93;
REAL8 tmp95=1.+tmp92+tmp94;
REAL8 tmp30=tmp1+tmp10+tmp12+tmp2+tmp3+tmp8;
REAL8 tmp14=sqrt(tmp13);
REAL8 tmp100=tmp30*tmp30;
REAL8 tmp101=-(tmp13*tmp4*tmp53*tmp60*tmp95);
REAL8 tmp102=tmp100+tmp101;
REAL8 tmp17=-(tmp15*tmp16*tmp9*x->data[0]);
REAL8 tmp18=tmp15*tmp16*tmp7*x->data[1];
REAL8 tmp19=tmp17+tmp18;
REAL8 tmp20=p->data[2]*tmp19;
REAL8 tmp21=tmp11*tmp15*tmp16*x->data[0];
REAL8 tmp22=-(tmp15*tmp16*tmp7*x->data[2]);
REAL8 tmp23=tmp21+tmp22;
REAL8 tmp24=p->data[1]*tmp23;
REAL8 tmp25=-(tmp11*tmp15*tmp16*x->data[1]);
REAL8 tmp26=tmp15*tmp16*tmp9*x->data[2];
REAL8 tmp27=tmp25+tmp26;
REAL8 tmp28=p->data[0]*tmp27;
REAL8 tmp29=tmp20+tmp24+tmp28;
REAL8 tmp120=1/tmp102;
REAL8 tmp124=tmp13*tmp58;
REAL8 tmp125=tmp1+tmp124+tmp2+tmp3;
REAL8 tmp126=1/tmp125;
REAL8 tmp105=tmp10*tmp16*tmp86*x->data[0];
REAL8 tmp106=-(tmp15*tmp16*x->data[0]);
REAL8 tmp107=-(tmp16*tmp7*tmp86*tmp9*x->data[1]);
REAL8 tmp108=tmp105+tmp106+tmp107;
REAL8 tmp110=-(tmp11*tmp16*tmp86*tmp9*x->data[0]);
REAL8 tmp111=tmp16*tmp7*tmp86*tmp9*x->data[2];
REAL8 tmp112=tmp110+tmp111;
REAL8 tmp114=tmp11*tmp16*tmp86*tmp9*x->data[1];
REAL8 tmp115=-(tmp10*tmp16*tmp86*x->data[2]);
REAL8 tmp116=tmp15*tmp16*x->data[2];
REAL8 tmp117=tmp114+tmp115+tmp116;
REAL8 tmp130=1/tmp60;
REAL8 tmp87=-(tmp16*tmp7*tmp86*tmp9*x->data[0]);
REAL8 tmp88=-(tmp10*tmp16*tmp86*x->data[1]);
REAL8 tmp89=tmp15*tmp16*x->data[1];
REAL8 tmp90=-(tmp11*tmp16*tmp86*tmp9*x->data[2]);
REAL8 tmp91=tmp87+tmp88+tmp89+tmp90;
REAL8 tmp144=tmp16*tmp19*x->data[1];
REAL8 tmp145=-(tmp16*tmp23*x->data[2]);
REAL8 tmp146=tmp144+tmp145;
REAL8 tmp147=p->data[0]*tmp146;
REAL8 tmp148=tmp16*tmp23*x->data[0];
REAL8 tmp149=-(tmp16*tmp27*x->data[1]);
REAL8 tmp150=tmp148+tmp149;
REAL8 tmp151=p->data[2]*tmp150;
REAL8 tmp152=-(tmp16*tmp19*x->data[0]);
REAL8 tmp153=tmp16*tmp27*x->data[2];
REAL8 tmp154=tmp152+tmp153;
REAL8 tmp155=p->data[1]*tmp154;
REAL8 tmp156=tmp147+tmp151+tmp155;
REAL8 tmp162=tmp156*tmp156;
REAL8 tmp32=2.*c1k5*tmp9;
REAL8 tmp33=4.*c2k5*tmp13*tmp9;
REAL8 tmp34=tmp32+tmp33;
REAL8 tmp36=1.*tmp34*tmp35;
REAL8 tmp37=2.*c1k4*tmp9;
REAL8 tmp38=4.*c2k4*tmp13*tmp9;
REAL8 tmp39=tmp37+tmp38;
REAL8 tmp40=1.*tmp39*tmp5;
REAL8 tmp42=2.*c1k3*tmp41*tmp9;
REAL8 tmp44=2.*c1k2*tmp43*tmp9;
REAL8 tmp45=tmp36+tmp40+tmp42+tmp44;
REAL8 tmp84=1/tmp83;
REAL8 tmp166=p->data[0]*tmp16*x->data[0];
REAL8 tmp167=p->data[1]*tmp16*x->data[1];
REAL8 tmp168=p->data[2]*tmp16*x->data[2];
REAL8 tmp169=tmp166+tmp167+tmp168;
REAL8 tmp170=tmp169*tmp169;
REAL8 tmp158=2.*tmp13*tmp57*tmp91;
REAL8 tmp159=2.*tmp58*tmp9;
REAL8 tmp160=tmp158+tmp159;
REAL8 tmp161=(1.0/(tmp125*tmp125));
REAL8 tmp171=-3.*eta;
REAL8 tmp172=26.+tmp171;
REAL8 tmp173=2.*eta*tmp172*tmp41;
REAL8 tmp174=6.*eta*tmp43;
REAL8 tmp175=1.+tmp173+tmp174;
REAL8 tmp176=log(tmp175);
REAL8 tmp177=1.+tmp176;
REAL8 tmp181=4.+tmp171;
REAL8 tmp184=1./(tmp30*tmp30*tmp30*tmp30);
REAL8 tmp186=(tmp169*tmp169*tmp169*tmp169);
REAL8 tmp187=tmp177*tmp177;
REAL8 tmp183=(tmp4*tmp4*tmp4);
REAL8 tmp185=(tmp53*tmp53*tmp53*tmp53);
REAL8 tmp192=(tmp95*tmp95*tmp95*tmp95);
REAL8 tmp31=4.*tmp30*tmp9;
REAL8 tmp85=-(eta*tmp13*tmp4*tmp45*tmp53*tmp60*tmp84);
REAL8 tmp96=2.*tmp13*tmp4*tmp53*tmp57*tmp91*tmp95;
REAL8 tmp97=-2.*tmp13*tmp60*tmp9*tmp95;
REAL8 tmp98=-2.*tmp4*tmp53*tmp60*tmp9*tmp95;
REAL8 tmp99=tmp31+tmp85+tmp96+tmp97+tmp98;
REAL8 tmp103=(1.0/(tmp102*tmp102));
REAL8 tmp196=tmp29*tmp29;
REAL8 tmp164=(1.0/(tmp60*tmp60));
REAL8 tmp109=p->data[2]*tmp108;
REAL8 tmp113=p->data[1]*tmp112;
REAL8 tmp118=p->data[0]*tmp117;
REAL8 tmp119=tmp109+tmp113+tmp118;
REAL8 tmp123=1/tmp53;
REAL8 tmp127=1/tmp95;
REAL8 tmp128=tmp102*tmp123*tmp126*tmp127*tmp43;
REAL8 tmp202=tmp126*tmp130*tmp162*tmp4;
REAL8 tmp203=tmp126*tmp170*tmp177*tmp4*tmp53*tmp95;
REAL8 tmp204=2.*eta*tmp181*tmp183*tmp184*tmp185*tmp186*tmp187*tmp192;
REAL8 tmp205=tmp120*tmp125*tmp130*tmp196*tmp4;
REAL8 tmp206=1.+tmp202+tmp203+tmp204+tmp205;
REAL8 tmp220=sqrt(tmp4*tmp4*tmp4);
REAL8 tmp222=1/mass1;
REAL8 tmp223=mass2*s1Vec->data[0]*tmp222;
REAL8 tmp224=1/mass2;
REAL8 tmp225=mass1*s2Vec->data[0]*tmp224;
REAL8 tmp226=tmp223+tmp225;
REAL8 tmp219=sqrt(tmp4);
REAL8 tmp131=tmp108*tmp16*x->data[1];
REAL8 tmp132=-(tmp112*tmp16*x->data[2]);
REAL8 tmp133=tmp131+tmp132;
REAL8 tmp134=p->data[0]*tmp133;
REAL8 tmp135=tmp112*tmp16*x->data[0];
REAL8 tmp136=-(tmp117*tmp16*x->data[1]);
REAL8 tmp137=tmp135+tmp136;
REAL8 tmp138=p->data[2]*tmp137;
REAL8 tmp139=-(tmp108*tmp16*x->data[0]);
REAL8 tmp140=tmp117*tmp16*x->data[2];
REAL8 tmp141=tmp139+tmp140;
REAL8 tmp142=p->data[1]*tmp141;
REAL8 tmp143=tmp134+tmp138+tmp142;
REAL8 tmp157=2.*tmp126*tmp130*tmp143*tmp156*tmp4;
REAL8 tmp163=-(tmp130*tmp160*tmp161*tmp162*tmp4);
REAL8 tmp165=2.*tmp126*tmp162*tmp164*tmp4*tmp57*tmp91;
REAL8 tmp178=eta*tmp126*tmp170*tmp177*tmp4*tmp45*tmp53*tmp84;
REAL8 tmp179=-(tmp160*tmp161*tmp170*tmp177*tmp4*tmp53*tmp95);
REAL8 tmp180=2.*tmp126*tmp170*tmp177*tmp9*tmp95;
REAL8 tmp197=-(tmp103*tmp125*tmp130*tmp196*tmp4*tmp99);
REAL8 tmp198=tmp120*tmp130*tmp160*tmp196*tmp4;
REAL8 tmp199=2.*tmp120*tmp125*tmp164*tmp196*tmp4*tmp57*tmp91;
REAL8 tmp200=2.*tmp119*tmp120*tmp125*tmp130*tmp29*tmp4;
REAL8 tmp232=tmp157+tmp163+tmp165+tmp178+tmp179+tmp180+tmp197+tmp198+tmp199+tmp200;
REAL8 tmp182=eta*eta;
REAL8 tmp238=tmp53*tmp53;
REAL8 tmp190=tmp4*tmp4;
REAL8 tmp241=tmp95*tmp95;
REAL8 tmp248=-16.*eta;
REAL8 tmp249=21.*tmp182;
REAL8 tmp250=tmp248+tmp249;
REAL8 tmp254=0.+tmp202+tmp203+tmp205;
REAL8 tmp256=-47.*eta;
REAL8 tmp257=54.*tmp182;
REAL8 tmp258=tmp219*tmp250*tmp254;
REAL8 tmp259=tmp256+tmp257+tmp258;
REAL8 tmp237=(eta*eta*eta);
REAL8 tmp240=1./(tmp125*tmp125*tmp125);
REAL8 tmp272=-6.*eta;
REAL8 tmp273=39.*tmp182;
REAL8 tmp274=tmp272+tmp273;
REAL8 tmp277=16.*eta;
REAL8 tmp278=147.*tmp182;
REAL8 tmp279=tmp219*tmp254*tmp274;
REAL8 tmp280=tmp277+tmp278+tmp279;
REAL8 tmp289=mass2*s1Vec->data[2]*tmp222;
REAL8 tmp290=mass1*s2Vec->data[2]*tmp224;
REAL8 tmp291=tmp289+tmp290;
REAL8 tmp239=-720.*tmp161*tmp183*tmp186*tmp187*tmp237*tmp238*tmp45*tmp84*tmp95;
REAL8 tmp242=720.*tmp160*tmp182*tmp183*tmp186*tmp187*tmp238*tmp240*tmp241;
REAL8 tmp243=-1440.*tmp161*tmp182*tmp186*tmp187*tmp190*tmp241*tmp53*tmp9;
REAL8 tmp244=103.*eta;
REAL8 tmp245=-60.*tmp182;
REAL8 tmp246=tmp244+tmp245;
REAL8 tmp247=2.*tmp219*tmp232*tmp246;
REAL8 tmp251=6.*tmp126*tmp170*tmp177*tmp190*tmp232*tmp250*tmp53*tmp95;
REAL8 tmp252=3.*eta;
REAL8 tmp253=23.+tmp252;
REAL8 tmp255=2.*eta*tmp232*tmp253*tmp254*tmp4;
REAL8 tmp260=6.*eta*tmp126*tmp170*tmp177*tmp220*tmp259*tmp45*tmp53*tmp84;
REAL8 tmp261=-6.*tmp160*tmp161*tmp170*tmp177*tmp220*tmp259*tmp53*tmp95;
REAL8 tmp262=12.*tmp126*tmp170*tmp177*tmp219*tmp259*tmp9*tmp95;
REAL8 tmp263=tmp239+tmp242+tmp243+tmp247+tmp251+tmp255+tmp260+tmp261+tmp262;
REAL8 tmp265=1620.*tmp161*tmp183*tmp186*tmp187*tmp237*tmp238*tmp45*tmp84*tmp95;
REAL8 tmp266=-1620.*tmp160*tmp182*tmp183*tmp186*tmp187*tmp238*tmp240*tmp241;
REAL8 tmp267=3240.*tmp161*tmp182*tmp186*tmp187*tmp190*tmp241*tmp53*tmp9;
REAL8 tmp268=-109.*eta;
REAL8 tmp269=51.*tmp182;
REAL8 tmp270=tmp268+tmp269;
REAL8 tmp271=4.*tmp219*tmp232*tmp270;
REAL8 tmp275=-6.*tmp126*tmp170*tmp177*tmp190*tmp232*tmp274*tmp53*tmp95;
REAL8 tmp276=-90.*eta*tmp232*tmp254*tmp4;
REAL8 tmp281=-6.*eta*tmp126*tmp170*tmp177*tmp220*tmp280*tmp45*tmp53*tmp84;
REAL8 tmp282=6.*tmp160*tmp161*tmp170*tmp177*tmp220*tmp280*tmp53*tmp95;
REAL8 tmp283=-12.*tmp126*tmp170*tmp177*tmp219*tmp280*tmp9*tmp95;
REAL8 tmp284=tmp265+tmp266+tmp267+tmp271+tmp275+tmp276+tmp281+tmp282+tmp283;
REAL8 tmp309=mass2*s1Vec->data[1]*tmp222;
REAL8 tmp310=mass1*s2Vec->data[1]*tmp224;
REAL8 tmp311=tmp309+tmp310;
REAL8 tmp331=tmp254*tmp254;
REAL8 tmp326=27.*eta;
REAL8 tmp327=-353.+tmp326;
REAL8 tmp328=2.*eta*tmp327;
REAL8 tmp329=-360.*tmp161*tmp182*tmp183*tmp186*tmp187*tmp238*tmp241;
REAL8 tmp330=2.*tmp219*tmp246*tmp254;
REAL8 tmp332=eta*tmp253*tmp331*tmp4;
REAL8 tmp333=6.*tmp126*tmp170*tmp177*tmp220*tmp259*tmp53*tmp95;
REAL8 tmp334=tmp328+tmp329+tmp330+tmp332+tmp333;
REAL8 tmp337=8.+tmp252;
REAL8 tmp338=-112.*eta*tmp337;
REAL8 tmp339=810.*tmp161*tmp182*tmp183*tmp186*tmp187*tmp238*tmp241;
REAL8 tmp340=4.*tmp219*tmp254*tmp270;
REAL8 tmp341=-45.*eta*tmp331*tmp4;
REAL8 tmp342=-6.*tmp126*tmp170*tmp177*tmp220*tmp280*tmp53*tmp95;
REAL8 tmp343=tmp338+tmp339+tmp340+tmp341+tmp342;
REAL8 tmp360=coeffs->d1v2*eta*tmp41*tmp9;
REAL8 tmp361=-8.*tmp9;
REAL8 tmp362=14.*tmp311;
REAL8 tmp363=-36.*tmp126*tmp170*tmp177*tmp220*tmp53*tmp9*tmp95;
REAL8 tmp364=-30.*tmp126*tmp170*tmp177*tmp220*tmp311*tmp53*tmp95;
REAL8 tmp365=3.*tmp219*tmp254*tmp9;
REAL8 tmp366=4.*tmp219*tmp254*tmp311;
REAL8 tmp367=tmp361+tmp362+tmp363+tmp364+tmp365+tmp366;
REAL8 tmp368=0.08333333333333333*eta*tmp16*tmp367;
REAL8 tmp369=-0.013888888888888888*tmp311*tmp334*tmp43;
REAL8 tmp370=0.006944444444444444*tmp343*tmp43*tmp9;
REAL8 tmp371=tmp309+tmp310+tmp360+tmp368+tmp369+tmp370;
REAL8 tmp389=sqrt(tmp125);
REAL8 tmp221=-36.*eta*tmp126*tmp170*tmp177*tmp220*tmp45*tmp53*tmp7*tmp84;
REAL8 tmp227=-30.*eta*tmp126*tmp170*tmp177*tmp220*tmp226*tmp45*tmp53*tmp84;
REAL8 tmp228=36.*tmp160*tmp161*tmp170*tmp177*tmp220*tmp53*tmp7*tmp95;
REAL8 tmp229=30.*tmp160*tmp161*tmp170*tmp177*tmp220*tmp226*tmp53*tmp95;
REAL8 tmp230=-72.*tmp126*tmp170*tmp177*tmp219*tmp7*tmp9*tmp95;
REAL8 tmp231=-60.*tmp126*tmp170*tmp177*tmp219*tmp226*tmp9*tmp95;
REAL8 tmp233=3.*tmp219*tmp232*tmp7;
REAL8 tmp234=4.*tmp219*tmp226*tmp232;
REAL8 tmp235=tmp221+tmp227+tmp228+tmp229+tmp230+tmp231+tmp233+tmp234;
REAL8 tmp236=0.08333333333333333*eta*tmp16*tmp235;
REAL8 tmp264=-0.013888888888888888*tmp226*tmp263*tmp43;
REAL8 tmp285=0.006944444444444444*tmp284*tmp43*tmp7;
REAL8 tmp286=tmp236+tmp264+tmp285;
REAL8 tmp287=tmp15*tmp286*tmp7;
REAL8 tmp288=-36.*eta*tmp11*tmp126*tmp170*tmp177*tmp220*tmp45*tmp53*tmp84;
REAL8 tmp292=-30.*eta*tmp126*tmp170*tmp177*tmp220*tmp291*tmp45*tmp53*tmp84;
REAL8 tmp293=36.*tmp11*tmp160*tmp161*tmp170*tmp177*tmp220*tmp53*tmp95;
REAL8 tmp294=30.*tmp160*tmp161*tmp170*tmp177*tmp220*tmp291*tmp53*tmp95;
REAL8 tmp295=-72.*tmp11*tmp126*tmp170*tmp177*tmp219*tmp9*tmp95;
REAL8 tmp296=-60.*tmp126*tmp170*tmp177*tmp219*tmp291*tmp9*tmp95;
REAL8 tmp297=3.*tmp11*tmp219*tmp232;
REAL8 tmp298=4.*tmp219*tmp232*tmp291;
REAL8 tmp299=tmp288+tmp292+tmp293+tmp294+tmp295+tmp296+tmp297+tmp298;
REAL8 tmp300=0.08333333333333333*eta*tmp16*tmp299;
REAL8 tmp301=-0.013888888888888888*tmp263*tmp291*tmp43;
REAL8 tmp302=0.006944444444444444*tmp11*tmp284*tmp43;
REAL8 tmp303=tmp300+tmp301+tmp302;
REAL8 tmp304=tmp11*tmp15*tmp303;
REAL8 tmp305=mass2*tmp222;
REAL8 tmp306=coeffs->d1v2*eta*tmp41;
REAL8 tmp307=14.*mass2*tmp222;
REAL8 tmp308=-36.*eta*tmp126*tmp170*tmp177*tmp220*tmp45*tmp53*tmp84*tmp9;
REAL8 tmp312=-30.*eta*tmp126*tmp170*tmp177*tmp220*tmp311*tmp45*tmp53*tmp84;
REAL8 tmp313=36.*tmp160*tmp161*tmp170*tmp177*tmp220*tmp53*tmp9*tmp95;
REAL8 tmp314=30.*tmp160*tmp161*tmp170*tmp177*tmp220*tmp311*tmp53*tmp95;
REAL8 tmp315=-72.*tmp10*tmp126*tmp170*tmp177*tmp219*tmp95;
REAL8 tmp316=-60.*tmp126*tmp170*tmp177*tmp219*tmp311*tmp9*tmp95;
REAL8 tmp317=-36.*tmp126*tmp170*tmp177*tmp220*tmp53*tmp95;
REAL8 tmp318=-30.*mass2*tmp126*tmp170*tmp177*tmp220*tmp222*tmp53*tmp95;
REAL8 tmp319=3.*tmp219*tmp232*tmp9;
REAL8 tmp320=4.*tmp219*tmp232*tmp311;
REAL8 tmp321=3.*tmp219*tmp254;
REAL8 tmp322=4.*mass2*tmp219*tmp222*tmp254;
REAL8 tmp323=-8.+tmp307+tmp308+tmp312+tmp313+tmp314+tmp315+tmp316+tmp317+tmp318+tmp319+tmp320+tmp321+tmp322;
REAL8 tmp324=0.08333333333333333*eta*tmp16*tmp323;
REAL8 tmp325=-0.013888888888888888*tmp263*tmp311*tmp43;
REAL8 tmp335=-0.013888888888888888*mass2*tmp222*tmp334*tmp43;
REAL8 tmp336=0.006944444444444444*tmp284*tmp43*tmp9;
REAL8 tmp344=0.006944444444444444*tmp343*tmp43;
REAL8 tmp345=tmp305+tmp306+tmp324+tmp325+tmp335+tmp336+tmp344;
REAL8 tmp346=tmp15*tmp345*tmp9;
REAL8 tmp347=coeffs->d1v2*eta*tmp41*tmp7;
REAL8 tmp348=-8.*tmp7;
REAL8 tmp349=14.*tmp226;
REAL8 tmp350=-36.*tmp126*tmp170*tmp177*tmp220*tmp53*tmp7*tmp95;
REAL8 tmp351=-30.*tmp126*tmp170*tmp177*tmp220*tmp226*tmp53*tmp95;
REAL8 tmp352=3.*tmp219*tmp254*tmp7;
REAL8 tmp353=4.*tmp219*tmp226*tmp254;
REAL8 tmp354=tmp348+tmp349+tmp350+tmp351+tmp352+tmp353;
REAL8 tmp355=0.08333333333333333*eta*tmp16*tmp354;
REAL8 tmp356=-0.013888888888888888*tmp226*tmp334*tmp43;
REAL8 tmp357=0.006944444444444444*tmp343*tmp43*tmp7;
REAL8 tmp358=tmp223+tmp225+tmp347+tmp355+tmp356+tmp357;
REAL8 tmp359=-(tmp358*tmp7*tmp86*tmp9);
REAL8 tmp372=-(tmp10*tmp371*tmp86);
REAL8 tmp373=tmp15*tmp371;
REAL8 tmp374=coeffs->d1v2*eta*tmp11*tmp41;
REAL8 tmp375=-8.*tmp11;
REAL8 tmp376=14.*tmp291;
REAL8 tmp377=-36.*tmp11*tmp126*tmp170*tmp177*tmp220*tmp53*tmp95;
REAL8 tmp378=-30.*tmp126*tmp170*tmp177*tmp220*tmp291*tmp53*tmp95;
REAL8 tmp379=3.*tmp11*tmp219*tmp254;
REAL8 tmp380=4.*tmp219*tmp254*tmp291;
REAL8 tmp381=tmp375+tmp376+tmp377+tmp378+tmp379+tmp380;
REAL8 tmp382=0.08333333333333333*eta*tmp16*tmp381;
REAL8 tmp383=-0.013888888888888888*tmp291*tmp334*tmp43;
REAL8 tmp384=0.006944444444444444*tmp11*tmp343*tmp43;
REAL8 tmp385=tmp289+tmp290+tmp374+tmp382+tmp383+tmp384;
REAL8 tmp386=-(tmp11*tmp385*tmp86*tmp9);
REAL8 tmp387=tmp287+tmp304+tmp346+tmp359+tmp372+tmp373+tmp386;
REAL8 tmp400=tmp15*tmp358*tmp7;
REAL8 tmp401=tmp15*tmp371*tmp9;
REAL8 tmp402=tmp11*tmp15*tmp385;
REAL8 tmp403=tmp400+tmp401+tmp402;
REAL8 tmp390=tmp4*tmp53*tmp95;
REAL8 tmp391=sqrt(tmp390);
REAL8 tmp392=-tmp391;
REAL8 tmp393=tmp120*tmp125*tmp4*tmp53*tmp95;
REAL8 tmp394=sqrt(tmp393);
REAL8 tmp395=tmp389*tmp394;
REAL8 tmp396=tmp392+tmp395;
REAL8 tmp397=1.+tmp202+tmp203+tmp205;
REAL8 tmp398=1./sqrt(tmp397);
REAL8 tmp409=1./sqrt(tmp125);
REAL8 tmp413=1./sqrt(tmp390);
REAL8 tmp419=1./sqrt(tmp393);
REAL8 tmp445=sqrt(tmp397);
REAL8 tmp435=tmp16*tmp358*x->data[0];
REAL8 tmp436=tmp16*tmp371*x->data[1];
REAL8 tmp437=tmp16*tmp385*x->data[2];
REAL8 tmp438=tmp435+tmp436+tmp437;
REAL8 tmp442=sqrt(tmp13*tmp13*tmp13);
REAL8 tmp443=(1.0/sqrt(tmp125*tmp125*tmp125));
REAL8 tmp406=(1.0/sqrt(tmp397*tmp397*tmp397));
REAL8 tmp446=1.+tmp445;
REAL8 tmp448=tmp125*tmp125;
REAL8 tmp449=-(tmp120*tmp190*tmp196*tmp448*tmp53*tmp95);
REAL8 tmp450=tmp162*tmp4;
REAL8 tmp451=1.+tmp202+tmp203+tmp205+tmp445;
REAL8 tmp452=-(tmp125*tmp451*tmp60);
REAL8 tmp453=tmp450+tmp452;
REAL8 tmp454=tmp4*tmp453*tmp53*tmp95;
REAL8 tmp455=tmp449+tmp454;
REAL8 tmp456=tmp438*tmp455;
REAL8 tmp457=tmp177*tmp4*tmp53*tmp95;
REAL8 tmp458=sqrt(tmp457);
REAL8 tmp459=tmp27*tmp358;
REAL8 tmp460=tmp23*tmp371;
REAL8 tmp461=tmp19*tmp385;
REAL8 tmp462=tmp459+tmp460+tmp461;
REAL8 tmp463=-(tmp219*tmp29*tmp389*tmp394*tmp462);
REAL8 tmp464=tmp146*tmp358;
REAL8 tmp465=tmp154*tmp371;
REAL8 tmp466=tmp150*tmp385;
REAL8 tmp467=tmp464+tmp465+tmp466;
REAL8 tmp468=tmp156*tmp219*tmp391*tmp467;
REAL8 tmp469=tmp463+tmp468;
REAL8 tmp470=-(tmp169*tmp391*tmp458*tmp469);
REAL8 tmp471=tmp456+tmp470;
REAL8 tmp473=1/tmp446;
REAL8 tmp414=eta*tmp4*tmp45*tmp53*tmp84;
REAL8 tmp415=2.*tmp9*tmp95;
REAL8 tmp416=tmp414+tmp415;
REAL8 tmp420=-(tmp103*tmp125*tmp4*tmp53*tmp95*tmp99);
REAL8 tmp421=eta*tmp120*tmp125*tmp4*tmp45*tmp53*tmp84;
REAL8 tmp422=tmp120*tmp160*tmp4*tmp53*tmp95;
REAL8 tmp423=2.*tmp120*tmp125*tmp9*tmp95;
REAL8 tmp424=tmp420+tmp421+tmp422+tmp423;
REAL8 tmp431=tmp16*tmp286*x->data[0];
REAL8 tmp432=tmp16*tmp303*x->data[2];
REAL8 tmp433=tmp16*tmp345*x->data[1];
REAL8 tmp434=tmp431+tmp432+tmp433;
REAL8 tmp538=coeffs->k5l*tmp81;
REAL8 tmp539=c0k5+tmp538+tmp61+tmp63;
REAL8 tmp511=tmp27*tmp286;
REAL8 tmp512=tmp19*tmp303;
REAL8 tmp513=tmp23*tmp345;
REAL8 tmp514=tmp117*tmp358;
REAL8 tmp515=tmp112*tmp371;
REAL8 tmp516=tmp108*tmp385;
REAL8 tmp517=tmp511+tmp512+tmp513+tmp514+tmp515+tmp516;
REAL8 tmp522=tmp146*tmp286;
REAL8 tmp523=tmp150*tmp303;
REAL8 tmp524=tmp154*tmp345;
REAL8 tmp525=tmp133*tmp358;
REAL8 tmp526=tmp141*tmp371;
REAL8 tmp527=tmp137*tmp385;
REAL8 tmp528=tmp522+tmp523+tmp524+tmp525+tmp526+tmp527;
REAL8 tmp506=1./sqrt(tmp457);
REAL8 tmp507=eta*tmp177*tmp4*tmp45*tmp53*tmp84;
REAL8 tmp508=2.*tmp177*tmp9*tmp95;
REAL8 tmp509=tmp507+tmp508;
REAL8 tmp494=0.5*tmp232*tmp398;
REAL8 tmp495=tmp157+tmp163+tmp165+tmp178+tmp179+tmp180+tmp197+tmp198+tmp199+tmp200+tmp494;
REAL8 tmp589=tmp156*tmp169*tmp219*tmp438*tmp458;
REAL8 tmp590=-(tmp170*tmp177*tmp4*tmp467*tmp53*tmp95);
REAL8 tmp591=tmp125*tmp451*tmp467;
REAL8 tmp592=tmp589+tmp590+tmp591;
REAL8 tmp444=1/tmp397;
REAL8 tmp536=2.*tmp102*tmp14;
REAL8 tmp537=4.*tmp219*tmp30;
REAL8 tmp540=1.*tmp35*tmp539;
REAL8 tmp541=1.+tmp540+tmp69+tmp72+tmp75+tmp79;
REAL8 tmp542=1/tmp541;
REAL8 tmp543=-2.*m1PlusEtaKK*tmp78;
REAL8 tmp544=2.*tmp74;
REAL8 tmp545=3.*tmp71;
REAL8 tmp546=4.*tmp68;
REAL8 tmp547=5.*tmp16*tmp539;
REAL8 tmp548=tmp546+tmp547;
REAL8 tmp549=1.*tmp16*tmp548;
REAL8 tmp550=tmp545+tmp549;
REAL8 tmp551=1.*tmp16*tmp550;
REAL8 tmp552=tmp544+tmp551;
REAL8 tmp553=1.*tmp16*tmp552;
REAL8 tmp554=tmp543+tmp553;
REAL8 tmp555=-(eta*tmp53*tmp542*tmp554);
REAL8 tmp556=2.*tmp219*tmp53*tmp95;
REAL8 tmp557=1.*tmp51;
REAL8 tmp558=1.*tmp13*tmp16;
REAL8 tmp559=tmp557+tmp558;
REAL8 tmp560=-2.*tmp559*tmp95;
REAL8 tmp561=tmp555+tmp556+tmp560;
REAL8 tmp562=-(tmp13*tmp561*tmp60);
REAL8 tmp563=tmp537+tmp562;
REAL8 tmp564=-2.*tmp14*tmp219*tmp563;
REAL8 tmp565=tmp536+tmp564;
REAL8 tmp447=(1.0/(tmp446*tmp446));
REAL8 tmp598=-(tmp156*tmp29*tmp389*tmp391*tmp394*tmp4*tmp462);
REAL8 tmp599=tmp120*tmp190*tmp196*tmp448*tmp467*tmp53*tmp95;
REAL8 tmp600=tmp4*tmp53*tmp592*tmp60*tmp95;
REAL8 tmp601=tmp598+tmp599+tmp600;
REAL8 tmp475=1./(tmp102*tmp102*tmp102);
REAL8 tmp478=pow(tmp125,-2.5);
REAL8 tmp483=(1.0/sqrt(tmp390*tmp390*tmp390));
REAL8 tmp485=(1.0/sqrt(tmp393*tmp393*tmp393));
REAL8 tmp649=tmp126*tmp219;
REAL8 tmp650=-tmp506;
REAL8 tmp651=tmp649+tmp650;
REAL8 tmp638=-(tmp4*tmp53*tmp95);
REAL8 tmp639=tmp1+tmp10+tmp12+tmp2+tmp3+tmp638+tmp8;
REAL8 tmp659=-4.*tmp220*tmp53*tmp95;
REAL8 tmp660=tmp30*tmp561;
REAL8 tmp661=tmp659+tmp660;
REAL8 tmp662=0.5*tmp120*tmp123*tmp127*tmp30*tmp43*tmp661;
REAL8 tmp663=tmp649+tmp662;
REAL8 tmp640=2.*tmp445;
REAL8 tmp641=1.+tmp640;
REAL8 tmp212=(1.0/(tmp95*tmp95));
REAL8 tmp642=-(tmp120*tmp13*tmp219*tmp29*tmp30*tmp391*tmp394*tmp409*tmp438*tmp57*tmp60*tmp639*tmp641);
REAL8 tmp643=tmp177*tmp190*tmp238*tmp241;
REAL8 tmp644=1./sqrt(tmp643);
REAL8 tmp645=-2.*tmp4*tmp53*tmp95;
REAL8 tmp646=tmp458*tmp561;
REAL8 tmp647=tmp645+tmp646;
REAL8 tmp648=-0.5*tmp219*tmp29*tmp389*tmp394*tmp446*tmp467*tmp644*tmp647;
REAL8 tmp652=tmp156*tmp219*tmp391*tmp462*tmp651;
REAL8 tmp653=-(tmp126*tmp13*tmp169*tmp57*tmp60);
REAL8 tmp654=tmp156*tmp219*tmp651;
REAL8 tmp655=-(tmp126*tmp13*tmp57);
REAL8 tmp656=tmp120*tmp126*tmp13*tmp30*tmp57*tmp639;
REAL8 tmp657=tmp655+tmp656;
REAL8 tmp658=tmp169*tmp60*tmp657;
REAL8 tmp664=-(tmp156*tmp219*tmp663);
REAL8 tmp665=tmp654+tmp658+tmp664;
REAL8 tmp666=tmp445*tmp665;
REAL8 tmp667=tmp653+tmp666;
REAL8 tmp668=tmp391*tmp462*tmp667;
REAL8 tmp669=tmp219*tmp29*tmp389*tmp394*tmp467*tmp641*tmp663;
REAL8 tmp670=tmp652+tmp668+tmp669;
REAL8 tmp671=tmp391*tmp670;
REAL8 tmp672=tmp648+tmp671;
REAL8 tmp673=tmp458*tmp672;
REAL8 tmp674=tmp642+tmp673;
REAL8 tmp676=1/tmp451;
REAL8 tmp215=(1.0/(tmp53*tmp53));
REAL8 tmp606=2.*eta*tmp219*tmp45*tmp53*tmp84;
REAL8 tmp607=-2.*eta*tmp45*tmp559*tmp84;
REAL8 tmp608=4.*c1k2*tmp9;
REAL8 tmp609=6.*c1k3*tmp9;
REAL8 tmp610=4.*tmp39;
REAL8 tmp611=5.*tmp16*tmp34;
REAL8 tmp612=tmp610+tmp611;
REAL8 tmp613=1.*tmp16*tmp612;
REAL8 tmp614=tmp609+tmp613;
REAL8 tmp615=1.*tmp16*tmp614;
REAL8 tmp616=tmp608+tmp615;
REAL8 tmp617=-(eta*tmp16*tmp53*tmp542*tmp616);
REAL8 tmp618=(1.0/(tmp541*tmp541));
REAL8 tmp619=eta*tmp45*tmp53*tmp554*tmp618;
REAL8 tmp620=-2.*eta*tmp43*tmp542*tmp554*tmp9;
REAL8 tmp621=0.+tmp606+tmp607+tmp617+tmp619+tmp620;
REAL8 tmp685=2.*tmp9;
REAL8 tmp686=-(eta*tmp4*tmp45*tmp53*tmp84);
REAL8 tmp687=-2.*tmp9*tmp95;
REAL8 tmp688=tmp685+tmp686+tmp687;
REAL8 tmp720=-(tmp160*tmp161*tmp219);
REAL8 tmp721=(1.0/sqrt(tmp457*tmp457*tmp457));
REAL8 tmp722=0.5*tmp509*tmp721;
REAL8 tmp723=tmp720+tmp722;
REAL8 tmp743=-4.*eta*tmp220*tmp45*tmp53*tmp84;
REAL8 tmp744=tmp30*tmp621;
REAL8 tmp745=-8.*tmp219*tmp9*tmp95;
REAL8 tmp746=2.*tmp561*tmp9;
REAL8 tmp747=tmp743+tmp744+tmp745+tmp746;
REAL8 tmp748=0.5*tmp120*tmp123*tmp127*tmp30*tmp43*tmp747;
REAL8 tmp749=-0.5*tmp103*tmp123*tmp127*tmp30*tmp43*tmp661*tmp99;
REAL8 tmp750=-0.5*eta*tmp120*tmp123*tmp212*tmp30*tmp43*tmp45*tmp661*tmp84;
REAL8 tmp751=-(tmp120*tmp127*tmp215*tmp30*tmp5*tmp661*tmp9);
REAL8 tmp752=1.*tmp120*tmp123*tmp127*tmp43*tmp661*tmp9;
REAL8 tmp753=tmp720+tmp748+tmp749+tmp750+tmp751+tmp752;
REAL8 tmp129=1./sqrt(tmp128);
REAL8 tmp210=sqrt(tmp206);
REAL8 ds010000=(1.*eta*(2.*tmp120*tmp14*tmp219*tmp387+tmp120*tmp130*tmp219*tmp29*tmp387*tmp389*tmp396*tmp398+2.*tmp119*tmp120*tmp14*tmp4+tmp119*tmp120*tmp130*tmp219*tmp389*tmp396*tmp398*tmp403-0.5*tmp120*tmp130*tmp219*tmp232*tmp29*tmp389*tmp396*tmp403*tmp406+0.5*tmp120*tmp130*tmp160*tmp219*tmp29*tmp396*tmp398*tmp403*tmp409+tmp120*tmp130*tmp219*tmp29*tmp389*tmp398*tmp403*(0.5*tmp160*tmp394*tmp409-0.5*tmp413*tmp416+0.5*tmp389*tmp419*tmp424)-0.5*tmp41*(2.*tmp286*tmp358+2.*tmp345*tmp371+2.*tmp303*tmp385-6.*tmp434*tmp438)+2.*coeffs->dheffSSv2*eta*s1Vec->data[1]*tmp5-0.25*tmp103*tmp130*tmp232*tmp413*tmp419*tmp443*tmp444*tmp447*tmp458*tmp565*tmp601-0.25*tmp103*tmp130*tmp232*tmp406*tmp413*tmp419*tmp443*tmp458*tmp473*tmp565*tmp601-0.75*tmp103*tmp130*tmp160*tmp398*tmp413*tmp419*tmp458*tmp473*tmp478*tmp565*tmp601-0.25*tmp103*tmp130*tmp398*tmp416*tmp419*tmp443*tmp458*tmp473*tmp483*tmp565*tmp601-0.25*tmp103*tmp130*tmp398*tmp413*tmp424*tmp443*tmp458*tmp473*tmp485*tmp565*tmp601+0.25*tmp103*tmp130*tmp398*tmp413*tmp419*tmp443*tmp473*tmp506*tmp509*tmp565*tmp601-tmp123*tmp127*tmp130*tmp160*tmp161*tmp394*tmp43*tmp674*tmp676+0.5*tmp123*tmp126*tmp127*tmp130*tmp419*tmp424*tmp43*tmp674*tmp676-2.*eta*tmp103*tmp220*tmp398*tmp413*tmp419*tmp442*tmp443*tmp45*tmp471*tmp473*tmp53*tmp57*tmp84-eta*tmp123*tmp126*tmp130*tmp212*tmp394*tmp43*tmp45*tmp674*tmp676*tmp84+2.*tmp120*tmp15*tmp29*tmp4*tmp9+2.*tmp120*tmp15*tmp219*tmp403*tmp9-2.*tmp126*tmp127*tmp130*tmp215*tmp394*tmp5*tmp674*tmp676*tmp9+2.*tmp120*tmp164*tmp219*tmp29*tmp389*tmp396*tmp398*tmp403*tmp57*tmp91+1.*tmp103*tmp164*tmp398*tmp413*tmp419*tmp443*tmp458*tmp473*tmp565*tmp57*tmp601*tmp91+2.*tmp123*tmp126*tmp127*tmp164*tmp394*tmp43*tmp57*tmp674*tmp676*tmp91+1.*tmp103*tmp220*tmp232*tmp413*tmp419*tmp442*tmp443*tmp444*tmp447*tmp471*tmp53*tmp57*tmp95+1.*tmp103*tmp220*tmp232*tmp406*tmp413*tmp419*tmp442*tmp443*tmp471*tmp473*tmp53*tmp57*tmp95+3.*tmp103*tmp160*tmp220*tmp398*tmp413*tmp419*tmp442*tmp471*tmp473*tmp478*tmp53*tmp57*tmp95+1.*tmp103*tmp220*tmp398*tmp416*tmp419*tmp442*tmp443*tmp471*tmp473*tmp483*tmp53*tmp57*tmp95+1.*tmp103*tmp220*tmp398*tmp413*tmp424*tmp442*tmp443*tmp471*tmp473*tmp485*tmp53*tmp57*tmp95-4.*tmp103*tmp219*tmp398*tmp413*tmp419*tmp442*tmp443*tmp471*tmp473*tmp57*tmp9*tmp95-6.*tmp103*tmp14*tmp220*tmp398*tmp413*tmp419*tmp443*tmp471*tmp473*tmp53*tmp57*tmp9*tmp95-2.*tmp103*tmp220*tmp398*tmp413*tmp419*tmp442*tmp443*tmp471*tmp473*tmp53*tmp91*tmp95-2.*tmp103*tmp14*tmp29*tmp4*tmp99-2.*tmp103*tmp14*tmp219*tmp403*tmp99-tmp103*tmp130*tmp219*tmp29*tmp389*tmp396*tmp398*tmp403*tmp99-tmp130*tmp398*tmp413*tmp419*tmp443*tmp458*tmp473*tmp475*tmp565*tmp601*tmp99+4.*tmp220*tmp398*tmp413*tmp419*tmp442*tmp443*tmp471*tmp473*tmp475*tmp53*tmp57*tmp95*tmp99+0.5*tmp103*tmp130*tmp398*tmp413*tmp419*tmp443*tmp458*tmp473*tmp601*(2.*tmp102*tmp15*tmp9-2.*tmp15*tmp219*tmp563*tmp9-2.*tmp14*tmp219*(-(tmp13*tmp60*tmp621)+8.*tmp219*tmp9-2.*tmp561*tmp60*tmp9+2.*tmp13*tmp561*tmp57*tmp91)+2.*tmp14*tmp99)+0.5*tmp103*tmp130*tmp398*tmp413*tmp419*tmp443*tmp458*tmp473*tmp565*(-(tmp119*tmp156*tmp389*tmp391*tmp394*tmp4*tmp462)-tmp143*tmp29*tmp389*tmp391*tmp394*tmp4*tmp462-0.5*tmp156*tmp160*tmp29*tmp391*tmp394*tmp4*tmp409*tmp462-0.5*tmp156*tmp29*tmp389*tmp394*tmp4*tmp413*tmp416*tmp462-0.5*tmp156*tmp29*tmp389*tmp391*tmp4*tmp419*tmp424*tmp462-tmp156*tmp29*tmp389*tmp391*tmp394*tmp4*tmp517+eta*tmp120*tmp190*tmp196*tmp448*tmp45*tmp467*tmp53*tmp84+eta*tmp4*tmp45*tmp53*tmp592*tmp60*tmp84+2.*tmp120*tmp125*tmp160*tmp190*tmp196*tmp467*tmp53*tmp95+2.*tmp119*tmp120*tmp190*tmp29*tmp448*tmp467*tmp53*tmp95+tmp120*tmp190*tmp196*tmp448*tmp528*tmp53*tmp95+2.*tmp120*tmp196*tmp4*tmp448*tmp467*tmp9*tmp95+2.*tmp592*tmp60*tmp9*tmp95-2.*tmp4*tmp53*tmp57*tmp592*tmp91*tmp95+tmp4*tmp53*tmp60*tmp95*(tmp156*tmp169*tmp219*tmp434*tmp458+tmp143*tmp169*tmp219*tmp438*tmp458+tmp160*tmp451*tmp467+tmp125*tmp467*tmp495+0.5*tmp156*tmp169*tmp219*tmp438*tmp506*tmp509+tmp125*tmp451*tmp528-eta*tmp170*tmp177*tmp4*tmp45*tmp467*tmp53*tmp84-tmp170*tmp177*tmp4*tmp528*tmp53*tmp95-2.*tmp170*tmp177*tmp467*tmp9*tmp95)-tmp103*tmp190*tmp196*tmp448*tmp467*tmp53*tmp95*tmp99)-2.*tmp103*tmp220*tmp398*tmp413*tmp419*tmp442*tmp443*tmp473*tmp53*tmp57*tmp95*(tmp434*tmp455-0.5*tmp169*tmp413*tmp416*tmp458*tmp469-0.5*tmp169*tmp391*tmp469*tmp506*tmp509-tmp169*tmp391*tmp458*(-(tmp119*tmp219*tmp389*tmp394*tmp462)-0.5*tmp160*tmp219*tmp29*tmp394*tmp409*tmp462-0.5*tmp219*tmp29*tmp389*tmp419*tmp424*tmp462+tmp143*tmp219*tmp391*tmp467+0.5*tmp156*tmp219*tmp413*tmp416*tmp467-tmp219*tmp29*tmp389*tmp394*tmp517+tmp156*tmp219*tmp391*tmp528)+tmp438*(-(eta*tmp120*tmp190*tmp196*tmp448*tmp45*tmp53*tmp84)+eta*tmp4*tmp45*tmp453*tmp53*tmp84-2.*tmp120*tmp125*tmp160*tmp190*tmp196*tmp53*tmp95-2.*tmp119*tmp120*tmp190*tmp29*tmp448*tmp53*tmp95-2.*tmp120*tmp196*tmp4*tmp448*tmp9*tmp95+2.*tmp453*tmp9*tmp95+tmp4*tmp53*(2.*tmp143*tmp156*tmp4-tmp160*tmp451*tmp60-tmp125*tmp495*tmp60+2.*tmp125*tmp451*tmp57*tmp91)*tmp95+tmp103*tmp190*tmp196*tmp448*tmp53*tmp95*tmp99))-0.5*tmp210*(-(tmp102*tmp123*tmp127*tmp160*tmp161*tmp43)-eta*tmp102*tmp123*tmp126*tmp212*tmp43*tmp45*tmp84-2.*tmp102*tmp126*tmp127*tmp215*tmp5*tmp9+tmp123*tmp126*tmp127*tmp43*tmp99)*(1.0/sqrt(tmp128*tmp128*tmp128))-tmp123*tmp126*tmp127*tmp130*tmp394*tmp43*tmp495*tmp674*(1.0/(tmp451*tmp451))+tmp123*tmp126*tmp127*tmp130*tmp394*tmp43*tmp676*(-(tmp120*tmp13*tmp219*tmp232*tmp29*tmp30*tmp391*tmp394*tmp398*tmp409*tmp438*tmp57*tmp60*tmp639)-tmp120*tmp13*tmp219*tmp29*tmp30*tmp391*tmp394*tmp409*tmp434*tmp57*tmp60*tmp639*tmp641-tmp119*tmp120*tmp13*tmp219*tmp30*tmp391*tmp394*tmp409*tmp438*tmp57*tmp60*tmp639*tmp641-0.5*tmp120*tmp13*tmp219*tmp29*tmp30*tmp394*tmp409*tmp413*tmp416*tmp438*tmp57*tmp60*tmp639*tmp641-0.5*tmp120*tmp13*tmp219*tmp29*tmp30*tmp391*tmp409*tmp419*tmp424*tmp438*tmp57*tmp60*tmp639*tmp641+0.5*tmp120*tmp13*tmp160*tmp219*tmp29*tmp30*tmp391*tmp394*tmp438*tmp443*tmp57*tmp60*tmp639*tmp641+0.5*tmp506*tmp509*tmp672-tmp120*tmp13*tmp219*tmp29*tmp30*tmp391*tmp394*tmp409*tmp438*tmp57*tmp60*tmp641*tmp688-2.*tmp120*tmp13*tmp219*tmp29*tmp391*tmp394*tmp409*tmp438*tmp57*tmp60*tmp639*tmp641*tmp9-2.*tmp120*tmp219*tmp29*tmp30*tmp391*tmp394*tmp409*tmp438*tmp57*tmp60*tmp639*tmp641*tmp9+2.*tmp120*tmp13*tmp219*tmp29*tmp30*tmp391*tmp394*tmp409*tmp438*tmp58*tmp639*tmp641*tmp91-tmp120*tmp13*tmp219*tmp29*tmp30*tmp391*tmp394*tmp409*tmp438*tmp60*tmp639*tmp641*tmp91+tmp103*tmp13*tmp219*tmp29*tmp30*tmp391*tmp394*tmp409*tmp438*tmp57*tmp60*tmp639*tmp641*tmp99+tmp458*(-0.25*tmp219*tmp232*tmp29*tmp389*tmp394*tmp398*tmp467*tmp644*tmp647-0.5*tmp119*tmp219*tmp389*tmp394*tmp446*tmp467*tmp644*tmp647-0.25*tmp160*tmp219*tmp29*tmp394*tmp409*tmp446*tmp467*tmp644*tmp647-0.25*tmp219*tmp29*tmp389*tmp419*tmp424*tmp446*tmp467*tmp644*tmp647-0.5*tmp219*tmp29*tmp389*tmp394*tmp446*tmp528*tmp644*tmp647+0.5*tmp413*tmp416*tmp670-0.5*tmp219*tmp29*tmp389*tmp394*tmp446*tmp467*tmp644*(0.5*tmp506*tmp509*tmp561+tmp458*tmp621-2.*eta*tmp4*tmp45*tmp53*tmp84-4.*tmp9*tmp95)+tmp391*(tmp143*tmp219*tmp391*tmp462*tmp651+0.5*tmp156*tmp219*tmp413*tmp416*tmp462*tmp651+tmp156*tmp219*tmp391*tmp517*tmp651+1.*tmp219*tmp232*tmp29*tmp389*tmp394*tmp398*tmp467*tmp663+tmp119*tmp219*tmp389*tmp394*tmp467*tmp641*tmp663+0.5*tmp160*tmp219*tmp29*tmp394*tmp409*tmp467*tmp641*tmp663+0.5*tmp219*tmp29*tmp389*tmp419*tmp424*tmp467*tmp641*tmp663+tmp219*tmp29*tmp389*tmp394*tmp528*tmp641*tmp663+0.5*tmp413*tmp416*tmp462*tmp667+tmp391*tmp517*tmp667+tmp156*tmp219*tmp391*tmp462*tmp723+tmp219*tmp29*tmp389*tmp394*tmp467*tmp641*tmp753+tmp391*tmp462*(tmp13*tmp160*tmp161*tmp169*tmp57*tmp60+0.5*tmp232*tmp398*tmp665-2.*tmp126*tmp169*tmp57*tmp60*tmp9+2.*tmp126*tmp13*tmp169*tmp58*tmp91-tmp126*tmp13*tmp169*tmp60*tmp91+tmp445*(tmp143*tmp219*tmp651-tmp143*tmp219*tmp663+tmp156*tmp219*tmp723-tmp156*tmp219*tmp753-2.*tmp169*tmp57*tmp657*tmp91+tmp169*tmp60*(tmp13*tmp160*tmp161*tmp57-tmp120*tmp13*tmp160*tmp161*tmp30*tmp57*tmp639+tmp120*tmp126*tmp13*tmp30*tmp57*tmp688-2.*tmp126*tmp57*tmp9+2.*tmp120*tmp126*tmp13*tmp57*tmp639*tmp9+2.*tmp120*tmp126*tmp30*tmp57*tmp639*tmp9-tmp126*tmp13*tmp91+tmp120*tmp126*tmp13*tmp30*tmp639*tmp91-tmp103*tmp126*tmp13*tmp30*tmp57*tmp639*tmp99))))+0.25*tmp219*tmp29*tmp389*tmp394*tmp446*tmp467*tmp647*(4.*tmp177*tmp241*tmp4*tmp53*tmp9+2.*eta*tmp177*tmp190*tmp238*tmp45*tmp84*tmp95)*(1.0/sqrt(tmp643*tmp643*tmp643))))+(0.5*tmp129*(tmp157+tmp163+tmp165+tmp178+tmp179+tmp180+tmp197+tmp198+tmp199+tmp200-16.*eta*tmp181*tmp183*tmp185*tmp186*tmp187*tmp192*tmp9*pow(tmp30,-5.)+16.*eta*tmp181*tmp184*tmp186*tmp187*tmp190*tmp192*tmp9*(tmp53*tmp53*tmp53)+8.*tmp181*tmp182*tmp183*tmp184*tmp185*tmp186*tmp187*tmp45*tmp84*(tmp95*tmp95*tmp95)))/sqrt(tmp206)))/sqrt(1.+2.*eta*(-1.+tmp129*tmp210+2.*tmp120*tmp14*tmp29*tmp4+2.*tmp120*tmp14*tmp219*tmp403+tmp120*tmp130*tmp219*tmp29*tmp389*tmp396*tmp398*tmp403+0.5*tmp103*tmp130*tmp398*tmp413*tmp419*tmp443*tmp458*tmp473*tmp565*tmp601+tmp123*tmp126*tmp127*tmp130*tmp394*tmp43*tmp674*tmp676-2.*tmp103*tmp220*tmp398*tmp413*tmp419*tmp442*tmp443*tmp471*tmp473*tmp53*tmp57*tmp95+coeffs->dheffSSv2*eta*tmp5*(s1Vec->data[0]*s1Vec->data[0]+s1Vec->data[1]*s1Vec->data[1]+s1Vec->data[2]*s1Vec->data[2]+s2Vec->data[0]*s2Vec->data[0]+s2Vec->data[1]*s2Vec->data[1]+s2Vec->data[2]*s2Vec->data[2])-0.5*tmp41*(tmp358*tmp358+tmp371*tmp371+tmp385*tmp385-3.*(tmp438*tmp438))));
