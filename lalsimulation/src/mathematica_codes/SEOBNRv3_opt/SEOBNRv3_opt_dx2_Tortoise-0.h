REAL8 tmp8=x->data[0]*x->data[0];
REAL8 tmp9=x->data[1]*x->data[1];
REAL8 tmp10=x->data[2]*x->data[2];
REAL8 tmp11=tmp10+tmp8+tmp9;
REAL8 tmp19=1./sqrt(tmp11);
REAL8 tmp14=sigmaKerr->data[0]*sigmaKerr->data[0];
REAL8 tmp15=sigmaKerr->data[1]*sigmaKerr->data[1];
REAL8 tmp16=sigmaKerr->data[2]*sigmaKerr->data[2];
REAL8 tmp17=tmp14+tmp15+tmp16;
REAL8 tmp35=coeffs->KK*eta;
REAL8 tmp36=-1.+tmp35;
REAL8 tmp51=pow(tmp11,-3.5);
REAL8 tmp12=1./(tmp11*tmp11*tmp11);
REAL8 tmp55=pow(tmp11,-2.5);
REAL8 tmp57=(1.0/(tmp11*tmp11));
REAL8 tmp59=(1.0/sqrt(tmp11*tmp11*tmp11));
REAL8 tmp39=1/tmp11;
REAL8 tmp61=1.*tmp19;
REAL8 tmp62=log(tmp61);
REAL8 tmp37=(1.0/(tmp36*tmp36));
REAL8 tmp38=1.*tmp37;
REAL8 tmp40=1.*tmp17*tmp39;
REAL8 tmp41=1/tmp36;
REAL8 tmp42=2.*tmp19*tmp41;
REAL8 tmp43=tmp38+tmp40+tmp42;
REAL8 tmp44=e3_x*tmp19*x->data[0];
REAL8 tmp45=e3_y*tmp19*x->data[1];
REAL8 tmp46=e3_z*tmp19*x->data[2];
REAL8 tmp47=tmp44+tmp45+tmp46;
REAL8 tmp65=1.*coeffs->k5*tmp55;
REAL8 tmp66=1.*coeffs->k4*tmp57;
REAL8 tmp67=1.*coeffs->k3*tmp59;
REAL8 tmp68=1.*coeffs->k2*tmp39;
REAL8 tmp69=1.*coeffs->k1*tmp19;
REAL8 tmp70=1.*coeffs->k5l*tmp55*tmp62;
REAL8 tmp71=1.+tmp65+tmp66+tmp67+tmp68+tmp69+tmp70;
REAL8 tmp48=tmp47*tmp47;
REAL8 tmp49=-tmp48;
REAL8 tmp50=1.+tmp49;
REAL8 tmp79=coeffs->k0*eta;
REAL8 tmp80=log(tmp71);
REAL8 tmp81=eta*tmp80;
REAL8 tmp82=1.+tmp79+tmp81;
REAL8 tmp33=tmp10+tmp14+tmp15+tmp16+tmp8+tmp9;
REAL8 tmp18=sqrt(tmp17);
REAL8 tmp90=tmp33*tmp33;
REAL8 tmp91=-(tmp11*tmp17*tmp43*tmp50*tmp82);
REAL8 tmp92=tmp90+tmp91;
REAL8 tmp20=-(e3_y*tmp19*x->data[0]);
REAL8 tmp21=e3_x*tmp19*x->data[1];
REAL8 tmp22=tmp20+tmp21;
REAL8 tmp23=p->data[2]*tmp22;
REAL8 tmp24=e3_z*tmp19*x->data[0];
REAL8 tmp25=-(e3_x*tmp19*x->data[2]);
REAL8 tmp26=tmp24+tmp25;
REAL8 tmp27=p->data[1]*tmp26;
REAL8 tmp28=-(e3_z*tmp19*x->data[1]);
REAL8 tmp29=e3_y*tmp19*x->data[2];
REAL8 tmp30=tmp28+tmp29;
REAL8 tmp31=p->data[0]*tmp30;
REAL8 tmp32=tmp23+tmp27+tmp31;
REAL8 tmp110=1/tmp92;
REAL8 tmp120=-3.*eta;
REAL8 tmp121=4.+tmp120;
REAL8 tmp127=p->data[0]*tmp19*x->data[0];
REAL8 tmp128=p->data[1]*tmp19*x->data[1];
REAL8 tmp129=p->data[2]*tmp19*x->data[2];
REAL8 tmp130=tmp127+tmp128+tmp129;
REAL8 tmp114=tmp17*tmp48;
REAL8 tmp115=tmp10+tmp114+tmp8+tmp9;
REAL8 tmp116=1/tmp115;
REAL8 tmp95=e3_y*tmp59*x->data[0]*x->data[2];
REAL8 tmp96=-(e3_x*tmp59*x->data[1]*x->data[2]);
REAL8 tmp97=tmp95+tmp96;
REAL8 tmp99=-(e3_z*tmp59*x->data[0]*x->data[2]);
REAL8 tmp100=e3_x*tmp10*tmp59;
REAL8 tmp101=-(e3_x*tmp19);
REAL8 tmp102=tmp100+tmp101+tmp99;
REAL8 tmp104=e3_z*tmp59*x->data[1]*x->data[2];
REAL8 tmp105=-(e3_y*tmp10*tmp59);
REAL8 tmp106=e3_y*tmp19;
REAL8 tmp107=tmp104+tmp105+tmp106;
REAL8 tmp74=-(e3_x*tmp59*x->data[0]*x->data[2]);
REAL8 tmp75=-(e3_y*tmp59*x->data[1]*x->data[2]);
REAL8 tmp76=-(e3_z*tmp10*tmp59);
REAL8 tmp77=e3_z*tmp19;
REAL8 tmp78=tmp74+tmp75+tmp76+tmp77;
REAL8 tmp135=1/tmp50;
REAL8 tmp157=tmp19*tmp22*x->data[1];
REAL8 tmp158=-(tmp19*tmp26*x->data[2]);
REAL8 tmp159=tmp157+tmp158;
REAL8 tmp160=p->data[0]*tmp159;
REAL8 tmp161=tmp19*tmp26*x->data[0];
REAL8 tmp162=-(tmp19*tmp30*x->data[1]);
REAL8 tmp163=tmp161+tmp162;
REAL8 tmp164=p->data[2]*tmp163;
REAL8 tmp165=-(tmp19*tmp22*x->data[0]);
REAL8 tmp166=tmp19*tmp30*x->data[2];
REAL8 tmp167=tmp165+tmp166;
REAL8 tmp168=p->data[1]*tmp167;
REAL8 tmp169=tmp160+tmp164+tmp168;
REAL8 tmp175=tmp169*tmp169;
REAL8 tmp52=-5.*coeffs->k5*tmp51*x->data[2];
REAL8 tmp53=-(coeffs->k5l*tmp51*x->data[2]);
REAL8 tmp54=-4.*coeffs->k4*tmp12*x->data[2];
REAL8 tmp56=-3.*coeffs->k3*tmp55*x->data[2];
REAL8 tmp58=-2.*coeffs->k2*tmp57*x->data[2];
REAL8 tmp60=-(coeffs->k1*tmp59*x->data[2]);
REAL8 tmp63=-5.*coeffs->k5l*tmp51*tmp62*x->data[2];
REAL8 tmp64=tmp52+tmp53+tmp54+tmp56+tmp58+tmp60+tmp63;
REAL8 tmp72=1/tmp71;
REAL8 tmp181=26.+tmp120;
REAL8 tmp182=2.*eta*tmp181*tmp59;
REAL8 tmp183=6.*eta*tmp39;
REAL8 tmp184=1.+tmp182+tmp183;
REAL8 tmp180=tmp130*tmp130;
REAL8 tmp171=2.*x->data[2];
REAL8 tmp172=2.*tmp17*tmp47*tmp78;
REAL8 tmp173=tmp171+tmp172;
REAL8 tmp174=(1.0/(tmp115*tmp115));
REAL8 tmp185=log(tmp184);
REAL8 tmp186=1.+tmp185;
REAL8 tmp122=-(p->data[0]*tmp59*x->data[0]*x->data[2]);
REAL8 tmp123=-(p->data[1]*tmp59*x->data[1]*x->data[2]);
REAL8 tmp124=-(p->data[2]*tmp10*tmp59);
REAL8 tmp125=p->data[2]*tmp19;
REAL8 tmp126=tmp122+tmp123+tmp124+tmp125;
REAL8 tmp84=-2.*tmp17*tmp57*x->data[2];
REAL8 tmp85=-2.*tmp41*tmp59*x->data[2];
REAL8 tmp86=tmp84+tmp85;
REAL8 tmp34=4.*tmp33*x->data[2];
REAL8 tmp73=-(eta*tmp11*tmp17*tmp43*tmp50*tmp64*tmp72);
REAL8 tmp83=2.*tmp11*tmp17*tmp43*tmp47*tmp78*tmp82;
REAL8 tmp87=-(tmp11*tmp17*tmp50*tmp82*tmp86);
REAL8 tmp88=-2.*tmp17*tmp43*tmp50*tmp82*x->data[2];
REAL8 tmp89=tmp34+tmp73+tmp83+tmp87+tmp88;
REAL8 tmp93=(1.0/(tmp92*tmp92));
REAL8 tmp197=tmp32*tmp32;
REAL8 tmp177=(1.0/(tmp50*tmp50));
REAL8 tmp98=p->data[2]*tmp97;
REAL8 tmp103=p->data[1]*tmp102;
REAL8 tmp108=p->data[0]*tmp107;
REAL8 tmp109=tmp103+tmp108+tmp98;
REAL8 tmp133=(tmp130*tmp130*tmp130*tmp130);
REAL8 tmp113=1/tmp43;
REAL8 tmp117=1/tmp82;
REAL8 tmp118=tmp113*tmp116*tmp117*tmp39*tmp92;
REAL8 tmp204=2.*eta*tmp121*tmp133*tmp39;
REAL8 tmp205=tmp11*tmp116*tmp135*tmp175;
REAL8 tmp206=tmp11*tmp116*tmp180*tmp186*tmp43*tmp82;
REAL8 tmp207=tmp11*tmp110*tmp115*tmp135*tmp197;
REAL8 tmp208=1.+tmp204+tmp205+tmp206+tmp207;
REAL8 tmp224=sqrt(tmp11*tmp11*tmp11);
REAL8 tmp188=-6.*eta*tmp181*tmp55*x->data[2];
REAL8 tmp189=-12.*eta*tmp57*x->data[2];
REAL8 tmp190=tmp188+tmp189;
REAL8 tmp191=1/tmp184;
REAL8 tmp222=sqrt(tmp11);
REAL8 tmp136=tmp19*tmp97*x->data[1];
REAL8 tmp137=-(tmp102*tmp19*x->data[2]);
REAL8 tmp138=-(tmp22*tmp59*x->data[1]*x->data[2]);
REAL8 tmp139=tmp10*tmp26*tmp59;
REAL8 tmp140=-(tmp19*tmp26);
REAL8 tmp141=tmp136+tmp137+tmp138+tmp139+tmp140;
REAL8 tmp142=p->data[0]*tmp141;
REAL8 tmp143=tmp102*tmp19*x->data[0];
REAL8 tmp144=-(tmp107*tmp19*x->data[1]);
REAL8 tmp145=-(tmp26*tmp59*x->data[0]*x->data[2]);
REAL8 tmp146=tmp30*tmp59*x->data[1]*x->data[2];
REAL8 tmp147=tmp143+tmp144+tmp145+tmp146;
REAL8 tmp148=p->data[2]*tmp147;
REAL8 tmp149=-(tmp19*tmp97*x->data[0]);
REAL8 tmp150=tmp107*tmp19*x->data[2];
REAL8 tmp151=tmp22*tmp59*x->data[0]*x->data[2];
REAL8 tmp152=-(tmp10*tmp30*tmp59);
REAL8 tmp153=tmp19*tmp30;
REAL8 tmp154=tmp149+tmp150+tmp151+tmp152+tmp153;
REAL8 tmp155=p->data[1]*tmp154;
REAL8 tmp156=tmp142+tmp148+tmp155;
REAL8 tmp170=2.*tmp11*tmp116*tmp135*tmp156*tmp169;
REAL8 tmp176=-(tmp11*tmp135*tmp173*tmp174*tmp175);
REAL8 tmp178=2.*tmp11*tmp116*tmp175*tmp177*tmp47*tmp78;
REAL8 tmp179=2.*tmp116*tmp135*tmp175*x->data[2];
REAL8 tmp187=eta*tmp11*tmp116*tmp180*tmp186*tmp43*tmp64*tmp72;
REAL8 tmp192=tmp11*tmp116*tmp180*tmp190*tmp191*tmp43*tmp82;
REAL8 tmp193=-(tmp11*tmp173*tmp174*tmp180*tmp186*tmp43*tmp82);
REAL8 tmp194=2.*tmp11*tmp116*tmp126*tmp130*tmp186*tmp43*tmp82;
REAL8 tmp195=tmp11*tmp116*tmp180*tmp186*tmp82*tmp86;
REAL8 tmp196=2.*tmp116*tmp180*tmp186*tmp43*tmp82*x->data[2];
REAL8 tmp198=-(tmp11*tmp115*tmp135*tmp197*tmp89*tmp93);
REAL8 tmp199=tmp11*tmp110*tmp135*tmp173*tmp197;
REAL8 tmp200=2.*tmp11*tmp110*tmp115*tmp177*tmp197*tmp47*tmp78;
REAL8 tmp201=2.*tmp109*tmp11*tmp110*tmp115*tmp135*tmp32;
REAL8 tmp202=2.*tmp110*tmp115*tmp135*tmp197*x->data[2];
REAL8 tmp237=tmp170+tmp176+tmp178+tmp179+tmp187+tmp192+tmp193+tmp194+tmp195+tmp196+tmp198+tmp199+tmp200+tmp201+tmp202;
REAL8 tmp240=0.+tmp205+tmp206+tmp207;
REAL8 tmp254=(tmp11*tmp11*tmp11);
REAL8 tmp255=tmp43*tmp43;
REAL8 tmp258=eta*eta;
REAL8 tmp256=tmp186*tmp186;
REAL8 tmp259=tmp82*tmp82;
REAL8 tmp131=(tmp130*tmp130*tmp130);
REAL8 tmp267=103.*eta;
REAL8 tmp268=-60.*tmp258;
REAL8 tmp269=tmp267+tmp268;
REAL8 tmp272=3.*eta;
REAL8 tmp273=23.+tmp272;
REAL8 tmp277=-16.*eta;
REAL8 tmp278=21.*tmp258;
REAL8 tmp279=tmp277+tmp278;
REAL8 tmp284=-47.*eta;
REAL8 tmp285=54.*tmp258;
REAL8 tmp286=tmp222*tmp240*tmp279;
REAL8 tmp287=tmp284+tmp285+tmp286;
REAL8 tmp275=tmp240*tmp240;
REAL8 tmp253=(eta*eta*eta);
REAL8 tmp261=1./(tmp115*tmp115*tmp115);
REAL8 tmp265=tmp11*tmp11;
REAL8 tmp311=-109.*eta;
REAL8 tmp312=51.*tmp258;
REAL8 tmp313=tmp311+tmp312;
REAL8 tmp318=-6.*eta;
REAL8 tmp319=39.*tmp258;
REAL8 tmp320=tmp318+tmp319;
REAL8 tmp325=16.*eta;
REAL8 tmp326=147.*tmp258;
REAL8 tmp327=tmp222*tmp240*tmp320;
REAL8 tmp328=tmp325+tmp326+tmp327;
REAL8 tmp257=-720.*tmp133*tmp174*tmp253*tmp254*tmp255*tmp256*tmp64*tmp72*tmp82;
REAL8 tmp260=-720.*tmp133*tmp174*tmp186*tmp190*tmp191*tmp254*tmp255*tmp258*tmp259;
REAL8 tmp262=720.*tmp133*tmp173*tmp254*tmp255*tmp256*tmp258*tmp259*tmp261;
REAL8 tmp263=-1440.*tmp126*tmp131*tmp174*tmp254*tmp255*tmp256*tmp258*tmp259;
REAL8 tmp264=-720.*tmp133*tmp174*tmp254*tmp256*tmp258*tmp259*tmp43*tmp86;
REAL8 tmp266=-2160.*tmp133*tmp174*tmp255*tmp256*tmp258*tmp259*tmp265*x->data[2];
REAL8 tmp270=2.*tmp222*tmp237*tmp269;
REAL8 tmp271=2.*tmp19*tmp240*tmp269*x->data[2];
REAL8 tmp274=2.*eta*tmp11*tmp237*tmp240*tmp273;
REAL8 tmp276=2.*eta*tmp273*tmp275*x->data[2];
REAL8 tmp280=tmp222*tmp237*tmp279;
REAL8 tmp281=tmp19*tmp240*tmp279*x->data[2];
REAL8 tmp282=tmp280+tmp281;
REAL8 tmp283=6.*tmp116*tmp180*tmp186*tmp224*tmp282*tmp43*tmp82;
REAL8 tmp288=6.*eta*tmp116*tmp180*tmp186*tmp224*tmp287*tmp43*tmp64*tmp72;
REAL8 tmp289=6.*tmp116*tmp180*tmp190*tmp191*tmp224*tmp287*tmp43*tmp82;
REAL8 tmp290=-6.*tmp173*tmp174*tmp180*tmp186*tmp224*tmp287*tmp43*tmp82;
REAL8 tmp291=12.*tmp116*tmp126*tmp130*tmp186*tmp224*tmp287*tmp43*tmp82;
REAL8 tmp292=6.*tmp116*tmp180*tmp186*tmp224*tmp287*tmp82*tmp86;
REAL8 tmp293=18.*tmp116*tmp180*tmp186*tmp222*tmp287*tmp43*tmp82*x->data[2];
REAL8 tmp294=tmp257+tmp260+tmp262+tmp263+tmp264+tmp266+tmp270+tmp271+tmp274+tmp276+tmp283+tmp288+tmp289+tmp290+tmp291+tmp292+tmp293;
REAL8 tmp296=27.*eta;
REAL8 tmp297=-353.+tmp296;
REAL8 tmp298=2.*eta*tmp297;
REAL8 tmp299=-360.*tmp133*tmp174*tmp254*tmp255*tmp256*tmp258*tmp259;
REAL8 tmp300=2.*tmp222*tmp240*tmp269;
REAL8 tmp301=eta*tmp11*tmp273*tmp275;
REAL8 tmp302=6.*tmp116*tmp180*tmp186*tmp224*tmp287*tmp43*tmp82;
REAL8 tmp303=tmp298+tmp299+tmp300+tmp301+tmp302;
REAL8 tmp305=1620.*tmp133*tmp174*tmp253*tmp254*tmp255*tmp256*tmp64*tmp72*tmp82;
REAL8 tmp306=1620.*tmp133*tmp174*tmp186*tmp190*tmp191*tmp254*tmp255*tmp258*tmp259;
REAL8 tmp307=-1620.*tmp133*tmp173*tmp254*tmp255*tmp256*tmp258*tmp259*tmp261;
REAL8 tmp308=3240.*tmp126*tmp131*tmp174*tmp254*tmp255*tmp256*tmp258*tmp259;
REAL8 tmp309=1620.*tmp133*tmp174*tmp254*tmp256*tmp258*tmp259*tmp43*tmp86;
REAL8 tmp310=4860.*tmp133*tmp174*tmp255*tmp256*tmp258*tmp259*tmp265*x->data[2];
REAL8 tmp314=4.*tmp222*tmp237*tmp313;
REAL8 tmp315=4.*tmp19*tmp240*tmp313*x->data[2];
REAL8 tmp316=-90.*eta*tmp11*tmp237*tmp240;
REAL8 tmp317=-90.*eta*tmp275*x->data[2];
REAL8 tmp321=tmp222*tmp237*tmp320;
REAL8 tmp322=tmp19*tmp240*tmp320*x->data[2];
REAL8 tmp323=tmp321+tmp322;
REAL8 tmp324=-6.*tmp116*tmp180*tmp186*tmp224*tmp323*tmp43*tmp82;
REAL8 tmp329=-6.*eta*tmp116*tmp180*tmp186*tmp224*tmp328*tmp43*tmp64*tmp72;
REAL8 tmp330=-6.*tmp116*tmp180*tmp190*tmp191*tmp224*tmp328*tmp43*tmp82;
REAL8 tmp331=6.*tmp173*tmp174*tmp180*tmp186*tmp224*tmp328*tmp43*tmp82;
REAL8 tmp332=-12.*tmp116*tmp126*tmp130*tmp186*tmp224*tmp328*tmp43*tmp82;
REAL8 tmp333=-6.*tmp116*tmp180*tmp186*tmp224*tmp328*tmp82*tmp86;
REAL8 tmp334=-18.*tmp116*tmp180*tmp186*tmp222*tmp328*tmp43*tmp82*x->data[2];
REAL8 tmp335=tmp305+tmp306+tmp307+tmp308+tmp309+tmp310+tmp314+tmp315+tmp316+tmp317+tmp324+tmp329+tmp330+tmp331+tmp332+tmp333+tmp334;
REAL8 tmp337=8.+tmp272;
REAL8 tmp338=-112.*eta*tmp337;
REAL8 tmp339=810.*tmp133*tmp174*tmp254*tmp255*tmp256*tmp258*tmp259;
REAL8 tmp340=4.*tmp222*tmp240*tmp313;
REAL8 tmp341=-45.*eta*tmp11*tmp275;
REAL8 tmp342=-6.*tmp116*tmp180*tmp186*tmp224*tmp328*tmp43*tmp82;
REAL8 tmp343=tmp338+tmp339+tmp340+tmp341+tmp342;
REAL8 tmp415=sqrt(tmp115);
REAL8 tmp223=-3.*coeffs->d1v2*eta*sigmaKerr->data[0]*tmp55*x->data[2];
REAL8 tmp225=-36.*eta*sigmaKerr->data[0]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp64*tmp72;
REAL8 tmp226=-30.*eta*sigmaStar->data[0]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp64*tmp72;
REAL8 tmp227=-36.*sigmaKerr->data[0]*tmp116*tmp180*tmp190*tmp191*tmp224*tmp43*tmp82;
REAL8 tmp228=-30.*sigmaStar->data[0]*tmp116*tmp180*tmp190*tmp191*tmp224*tmp43*tmp82;
REAL8 tmp229=36.*sigmaKerr->data[0]*tmp173*tmp174*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp230=30.*sigmaStar->data[0]*tmp173*tmp174*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp231=-72.*sigmaKerr->data[0]*tmp116*tmp126*tmp130*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp232=-60.*sigmaStar->data[0]*tmp116*tmp126*tmp130*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp233=-36.*sigmaKerr->data[0]*tmp116*tmp180*tmp186*tmp224*tmp82*tmp86;
REAL8 tmp234=-30.*sigmaStar->data[0]*tmp116*tmp180*tmp186*tmp224*tmp82*tmp86;
REAL8 tmp235=-108.*sigmaKerr->data[0]*tmp116*tmp180*tmp186*tmp222*tmp43*tmp82*x->data[2];
REAL8 tmp236=-90.*sigmaStar->data[0]*tmp116*tmp180*tmp186*tmp222*tmp43*tmp82*x->data[2];
REAL8 tmp238=3.*sigmaKerr->data[0]*tmp222*tmp237;
REAL8 tmp239=4.*sigmaStar->data[0]*tmp222*tmp237;
REAL8 tmp241=3.*sigmaKerr->data[0]*tmp19*tmp240*x->data[2];
REAL8 tmp242=4.*sigmaStar->data[0]*tmp19*tmp240*x->data[2];
REAL8 tmp243=tmp225+tmp226+tmp227+tmp228+tmp229+tmp230+tmp231+tmp232+tmp233+tmp234+tmp235+tmp236+tmp238+tmp239+tmp241+tmp242;
REAL8 tmp244=0.08333333333333333*eta*tmp19*tmp243;
REAL8 tmp245=-8.*sigmaKerr->data[0];
REAL8 tmp246=14.*sigmaStar->data[0];
REAL8 tmp247=-36.*sigmaKerr->data[0]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp248=-30.*sigmaStar->data[0]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp249=3.*sigmaKerr->data[0]*tmp222*tmp240;
REAL8 tmp250=4.*sigmaStar->data[0]*tmp222*tmp240;
REAL8 tmp251=tmp245+tmp246+tmp247+tmp248+tmp249+tmp250;
REAL8 tmp252=-0.08333333333333333*eta*tmp251*tmp59*x->data[2];
REAL8 tmp295=-0.013888888888888888*sigmaStar->data[0]*tmp294*tmp39;
REAL8 tmp304=0.027777777777777776*sigmaStar->data[0]*tmp303*tmp57*x->data[2];
REAL8 tmp336=0.006944444444444444*sigmaKerr->data[0]*tmp335*tmp39;
REAL8 tmp344=-0.013888888888888888*sigmaKerr->data[0]*tmp343*tmp57*x->data[2];
REAL8 tmp345=tmp223+tmp244+tmp252+tmp295+tmp304+tmp336+tmp344;
REAL8 tmp346=e3_x*tmp345;
REAL8 tmp347=-3.*coeffs->d1v2*eta*sigmaKerr->data[1]*tmp55*x->data[2];
REAL8 tmp348=-36.*eta*sigmaKerr->data[1]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp64*tmp72;
REAL8 tmp349=-30.*eta*sigmaStar->data[1]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp64*tmp72;
REAL8 tmp350=-36.*sigmaKerr->data[1]*tmp116*tmp180*tmp190*tmp191*tmp224*tmp43*tmp82;
REAL8 tmp351=-30.*sigmaStar->data[1]*tmp116*tmp180*tmp190*tmp191*tmp224*tmp43*tmp82;
REAL8 tmp352=36.*sigmaKerr->data[1]*tmp173*tmp174*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp353=30.*sigmaStar->data[1]*tmp173*tmp174*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp354=-72.*sigmaKerr->data[1]*tmp116*tmp126*tmp130*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp355=-60.*sigmaStar->data[1]*tmp116*tmp126*tmp130*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp356=-36.*sigmaKerr->data[1]*tmp116*tmp180*tmp186*tmp224*tmp82*tmp86;
REAL8 tmp357=-30.*sigmaStar->data[1]*tmp116*tmp180*tmp186*tmp224*tmp82*tmp86;
REAL8 tmp358=-108.*sigmaKerr->data[1]*tmp116*tmp180*tmp186*tmp222*tmp43*tmp82*x->data[2];
REAL8 tmp359=-90.*sigmaStar->data[1]*tmp116*tmp180*tmp186*tmp222*tmp43*tmp82*x->data[2];
REAL8 tmp360=3.*sigmaKerr->data[1]*tmp222*tmp237;
REAL8 tmp361=4.*sigmaStar->data[1]*tmp222*tmp237;
REAL8 tmp362=3.*sigmaKerr->data[1]*tmp19*tmp240*x->data[2];
REAL8 tmp363=4.*sigmaStar->data[1]*tmp19*tmp240*x->data[2];
REAL8 tmp364=tmp348+tmp349+tmp350+tmp351+tmp352+tmp353+tmp354+tmp355+tmp356+tmp357+tmp358+tmp359+tmp360+tmp361+tmp362+tmp363;
REAL8 tmp365=0.08333333333333333*eta*tmp19*tmp364;
REAL8 tmp366=-8.*sigmaKerr->data[1];
REAL8 tmp367=14.*sigmaStar->data[1];
REAL8 tmp368=-36.*sigmaKerr->data[1]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp369=-30.*sigmaStar->data[1]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp370=3.*sigmaKerr->data[1]*tmp222*tmp240;
REAL8 tmp371=4.*sigmaStar->data[1]*tmp222*tmp240;
REAL8 tmp372=tmp366+tmp367+tmp368+tmp369+tmp370+tmp371;
REAL8 tmp373=-0.08333333333333333*eta*tmp372*tmp59*x->data[2];
REAL8 tmp374=-0.013888888888888888*sigmaStar->data[1]*tmp294*tmp39;
REAL8 tmp375=0.027777777777777776*sigmaStar->data[1]*tmp303*tmp57*x->data[2];
REAL8 tmp376=0.006944444444444444*sigmaKerr->data[1]*tmp335*tmp39;
REAL8 tmp377=-0.013888888888888888*sigmaKerr->data[1]*tmp343*tmp57*x->data[2];
REAL8 tmp378=tmp347+tmp365+tmp373+tmp374+tmp375+tmp376+tmp377;
REAL8 tmp379=e3_y*tmp378;
REAL8 tmp380=-3.*coeffs->d1v2*eta*sigmaKerr->data[2]*tmp55*x->data[2];
REAL8 tmp381=-36.*eta*sigmaKerr->data[2]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp64*tmp72;
REAL8 tmp382=-30.*eta*sigmaStar->data[2]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp64*tmp72;
REAL8 tmp383=-36.*sigmaKerr->data[2]*tmp116*tmp180*tmp190*tmp191*tmp224*tmp43*tmp82;
REAL8 tmp384=-30.*sigmaStar->data[2]*tmp116*tmp180*tmp190*tmp191*tmp224*tmp43*tmp82;
REAL8 tmp385=36.*sigmaKerr->data[2]*tmp173*tmp174*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp386=30.*sigmaStar->data[2]*tmp173*tmp174*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp387=-72.*sigmaKerr->data[2]*tmp116*tmp126*tmp130*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp388=-60.*sigmaStar->data[2]*tmp116*tmp126*tmp130*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp389=-36.*sigmaKerr->data[2]*tmp116*tmp180*tmp186*tmp224*tmp82*tmp86;
REAL8 tmp390=-30.*sigmaStar->data[2]*tmp116*tmp180*tmp186*tmp224*tmp82*tmp86;
REAL8 tmp391=-108.*sigmaKerr->data[2]*tmp116*tmp180*tmp186*tmp222*tmp43*tmp82*x->data[2];
REAL8 tmp392=-90.*sigmaStar->data[2]*tmp116*tmp180*tmp186*tmp222*tmp43*tmp82*x->data[2];
REAL8 tmp393=3.*sigmaKerr->data[2]*tmp222*tmp237;
REAL8 tmp394=4.*sigmaStar->data[2]*tmp222*tmp237;
REAL8 tmp395=3.*sigmaKerr->data[2]*tmp19*tmp240*x->data[2];
REAL8 tmp396=4.*sigmaStar->data[2]*tmp19*tmp240*x->data[2];
REAL8 tmp397=tmp381+tmp382+tmp383+tmp384+tmp385+tmp386+tmp387+tmp388+tmp389+tmp390+tmp391+tmp392+tmp393+tmp394+tmp395+tmp396;
REAL8 tmp398=0.08333333333333333*eta*tmp19*tmp397;
REAL8 tmp399=-8.*sigmaKerr->data[2];
REAL8 tmp400=14.*sigmaStar->data[2];
REAL8 tmp401=-36.*sigmaKerr->data[2]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp402=-30.*sigmaStar->data[2]*tmp116*tmp180*tmp186*tmp224*tmp43*tmp82;
REAL8 tmp403=3.*sigmaKerr->data[2]*tmp222*tmp240;
REAL8 tmp404=4.*sigmaStar->data[2]*tmp222*tmp240;
REAL8 tmp405=tmp399+tmp400+tmp401+tmp402+tmp403+tmp404;
REAL8 tmp406=-0.08333333333333333*eta*tmp405*tmp59*x->data[2];
REAL8 tmp407=-0.013888888888888888*sigmaStar->data[2]*tmp294*tmp39;
REAL8 tmp408=0.027777777777777776*sigmaStar->data[2]*tmp303*tmp57*x->data[2];
REAL8 tmp409=0.006944444444444444*sigmaKerr->data[2]*tmp335*tmp39;
REAL8 tmp410=-0.013888888888888888*sigmaKerr->data[2]*tmp343*tmp57*x->data[2];
REAL8 tmp411=tmp380+tmp398+tmp406+tmp407+tmp408+tmp409+tmp410;
REAL8 tmp412=e3_z*tmp411;
REAL8 tmp413=tmp346+tmp379+tmp412;
REAL8 tmp426=coeffs->d1v2*eta*sigmaKerr->data[0]*tmp59;
REAL8 tmp427=0.08333333333333333*eta*tmp19*tmp251;
REAL8 tmp428=-0.013888888888888888*sigmaStar->data[0]*tmp303*tmp39;
REAL8 tmp429=0.006944444444444444*sigmaKerr->data[0]*tmp343*tmp39;
REAL8 tmp430=sigmaStar->data[0]+tmp426+tmp427+tmp428+tmp429;
REAL8 tmp431=e3_x*tmp430;
REAL8 tmp432=coeffs->d1v2*eta*sigmaKerr->data[1]*tmp59;
REAL8 tmp433=0.08333333333333333*eta*tmp19*tmp372;
REAL8 tmp434=-0.013888888888888888*sigmaStar->data[1]*tmp303*tmp39;
REAL8 tmp435=0.006944444444444444*sigmaKerr->data[1]*tmp343*tmp39;
REAL8 tmp436=sigmaStar->data[1]+tmp432+tmp433+tmp434+tmp435;
REAL8 tmp437=e3_y*tmp436;
REAL8 tmp438=coeffs->d1v2*eta*sigmaKerr->data[2]*tmp59;
REAL8 tmp439=0.08333333333333333*eta*tmp19*tmp405;
REAL8 tmp440=-0.013888888888888888*sigmaStar->data[2]*tmp303*tmp39;
REAL8 tmp441=0.006944444444444444*sigmaKerr->data[2]*tmp343*tmp39;
REAL8 tmp442=sigmaStar->data[2]+tmp438+tmp439+tmp440+tmp441;
REAL8 tmp443=e3_z*tmp442;
REAL8 tmp444=tmp431+tmp437+tmp443;
REAL8 tmp416=tmp11*tmp43*tmp82;
REAL8 tmp417=sqrt(tmp416);
REAL8 tmp418=-tmp417;
REAL8 tmp419=tmp11*tmp110*tmp115*tmp43*tmp82;
REAL8 tmp420=sqrt(tmp419);
REAL8 tmp421=tmp415*tmp420;
REAL8 tmp422=tmp418+tmp421;
REAL8 tmp423=1.+tmp205+tmp206+tmp207;
REAL8 tmp424=1./sqrt(tmp423);
REAL8 tmp450=1./sqrt(tmp115);
REAL8 tmp483=tmp19*tmp430*x->data[0];
REAL8 tmp484=tmp19*tmp436*x->data[1];
REAL8 tmp485=tmp19*tmp442*x->data[2];
REAL8 tmp486=tmp483+tmp484+tmp485;
REAL8 tmp455=1./sqrt(tmp416);
REAL8 tmp462=1./sqrt(tmp419);
REAL8 tmp500=sqrt(tmp423);
REAL8 tmp497=sqrt(tmp17*tmp17*tmp17);
REAL8 tmp498=(1.0/sqrt(tmp115*tmp115*tmp115));
REAL8 tmp447=(1.0/sqrt(tmp423*tmp423*tmp423));
REAL8 tmp501=1.+tmp500;
REAL8 tmp503=tmp115*tmp115;
REAL8 tmp504=-(tmp110*tmp197*tmp265*tmp43*tmp503*tmp82);
REAL8 tmp505=tmp11*tmp175;
REAL8 tmp506=1.+tmp205+tmp206+tmp207+tmp500;
REAL8 tmp507=-(tmp115*tmp50*tmp506);
REAL8 tmp508=tmp505+tmp507;
REAL8 tmp509=tmp11*tmp43*tmp508*tmp82;
REAL8 tmp510=tmp504+tmp509;
REAL8 tmp511=tmp486*tmp510;
REAL8 tmp512=tmp11*tmp186*tmp43*tmp82;
REAL8 tmp513=sqrt(tmp512);
REAL8 tmp514=tmp30*tmp430;
REAL8 tmp515=tmp26*tmp436;
REAL8 tmp516=tmp22*tmp442;
REAL8 tmp517=tmp514+tmp515+tmp516;
REAL8 tmp518=-(tmp222*tmp32*tmp415*tmp420*tmp517);
REAL8 tmp519=tmp159*tmp430;
REAL8 tmp520=tmp167*tmp436;
REAL8 tmp521=tmp163*tmp442;
REAL8 tmp522=tmp519+tmp520+tmp521;
REAL8 tmp523=tmp169*tmp222*tmp417*tmp522;
REAL8 tmp524=tmp518+tmp523;
REAL8 tmp525=-(tmp130*tmp417*tmp513*tmp524);
REAL8 tmp526=tmp511+tmp525;
REAL8 tmp528=1/tmp501;
REAL8 tmp456=eta*tmp11*tmp43*tmp64*tmp72;
REAL8 tmp457=tmp11*tmp82*tmp86;
REAL8 tmp458=2.*tmp43*tmp82*x->data[2];
REAL8 tmp459=tmp456+tmp457+tmp458;
REAL8 tmp463=-(tmp11*tmp115*tmp43*tmp82*tmp89*tmp93);
REAL8 tmp464=eta*tmp11*tmp110*tmp115*tmp43*tmp64*tmp72;
REAL8 tmp465=tmp11*tmp110*tmp173*tmp43*tmp82;
REAL8 tmp466=tmp11*tmp110*tmp115*tmp82*tmp86;
REAL8 tmp467=2.*tmp110*tmp115*tmp43*tmp82*x->data[2];
REAL8 tmp468=tmp463+tmp464+tmp465+tmp466+tmp467;
REAL8 tmp475=tmp19*tmp345*x->data[0];
REAL8 tmp476=tmp19*tmp378*x->data[1];
REAL8 tmp477=tmp19*tmp411*x->data[2];
REAL8 tmp478=-(tmp430*tmp59*x->data[0]*x->data[2]);
REAL8 tmp479=-(tmp436*tmp59*x->data[1]*x->data[2]);
REAL8 tmp480=-(tmp10*tmp442*tmp59);
REAL8 tmp481=tmp19*tmp442;
REAL8 tmp482=tmp475+tmp476+tmp477+tmp478+tmp479+tmp480+tmp481;
REAL8 tmp601=coeffs->k5l*tmp62;
REAL8 tmp602=coeffs->k5+tmp601;
REAL8 tmp572=tmp30*tmp345;
REAL8 tmp573=tmp26*tmp378;
REAL8 tmp574=tmp22*tmp411;
REAL8 tmp575=tmp107*tmp430;
REAL8 tmp576=tmp102*tmp436;
REAL8 tmp577=tmp442*tmp97;
REAL8 tmp578=tmp572+tmp573+tmp574+tmp575+tmp576+tmp577;
REAL8 tmp584=tmp159*tmp345;
REAL8 tmp585=tmp167*tmp378;
REAL8 tmp586=tmp163*tmp411;
REAL8 tmp587=tmp141*tmp430;
REAL8 tmp588=tmp154*tmp436;
REAL8 tmp589=tmp147*tmp442;
REAL8 tmp590=tmp584+tmp585+tmp586+tmp587+tmp588+tmp589;
REAL8 tmp565=1./sqrt(tmp512);
REAL8 tmp566=eta*tmp11*tmp186*tmp43*tmp64*tmp72;
REAL8 tmp567=tmp11*tmp190*tmp191*tmp43*tmp82;
REAL8 tmp568=tmp11*tmp186*tmp82*tmp86;
REAL8 tmp569=2.*tmp186*tmp43*tmp82*x->data[2];
REAL8 tmp570=tmp566+tmp567+tmp568+tmp569;
REAL8 tmp551=0.5*tmp237*tmp424;
REAL8 tmp552=tmp170+tmp176+tmp178+tmp179+tmp187+tmp192+tmp193+tmp194+tmp195+tmp196+tmp198+tmp199+tmp200+tmp201+tmp202+tmp551;
REAL8 tmp658=tmp130*tmp169*tmp222*tmp486*tmp513;
REAL8 tmp659=-(tmp11*tmp180*tmp186*tmp43*tmp522*tmp82);
REAL8 tmp660=tmp115*tmp506*tmp522;
REAL8 tmp661=tmp658+tmp659+tmp660;
REAL8 tmp499=1/tmp423;
REAL8 tmp599=2.*tmp18*tmp92;
REAL8 tmp600=4.*tmp222*tmp33;
REAL8 tmp603=1.*tmp55*tmp602;
REAL8 tmp604=1.+tmp603+tmp66+tmp67+tmp68+tmp69;
REAL8 tmp605=1/tmp604;
REAL8 tmp606=2.*coeffs->k2;
REAL8 tmp607=3.*coeffs->k3;
REAL8 tmp608=4.*coeffs->k4;
REAL8 tmp609=5.*tmp19*tmp602;
REAL8 tmp610=tmp608+tmp609;
REAL8 tmp611=1.*tmp19*tmp610;
REAL8 tmp612=tmp607+tmp611;
REAL8 tmp613=1.*tmp19*tmp612;
REAL8 tmp614=tmp606+tmp613;
REAL8 tmp615=1.*tmp19*tmp614;
REAL8 tmp616=coeffs->k1+tmp615;
REAL8 tmp617=-(eta*tmp43*tmp605*tmp616);
REAL8 tmp618=2.*tmp222*tmp43*tmp82;
REAL8 tmp619=1.*tmp41;
REAL8 tmp620=1.*tmp17*tmp19;
REAL8 tmp621=tmp619+tmp620;
REAL8 tmp622=-2.*tmp621*tmp82;
REAL8 tmp623=tmp617+tmp618+tmp622;
REAL8 tmp624=-(tmp17*tmp50*tmp623);
REAL8 tmp625=tmp600+tmp624;
REAL8 tmp626=-2.*tmp18*tmp222*tmp625;
REAL8 tmp627=tmp599+tmp626;
REAL8 tmp502=(1.0/(tmp501*tmp501));
REAL8 tmp668=-(tmp11*tmp169*tmp32*tmp415*tmp417*tmp420*tmp517);
REAL8 tmp669=tmp110*tmp197*tmp265*tmp43*tmp503*tmp522*tmp82;
REAL8 tmp670=tmp11*tmp43*tmp50*tmp661*tmp82;
REAL8 tmp671=tmp668+tmp669+tmp670;
REAL8 tmp530=1./(tmp92*tmp92*tmp92);
REAL8 tmp533=pow(tmp115,-2.5);
REAL8 tmp538=(1.0/sqrt(tmp416*tmp416*tmp416));
REAL8 tmp540=(1.0/sqrt(tmp419*tmp419*tmp419));
REAL8 tmp726=tmp116*tmp222;
REAL8 tmp727=-tmp565;
REAL8 tmp728=tmp726+tmp727;
REAL8 tmp715=-(tmp11*tmp43*tmp82);
REAL8 tmp716=tmp10+tmp14+tmp15+tmp16+tmp715+tmp8+tmp9;
REAL8 tmp736=-4.*tmp224*tmp43*tmp82;
REAL8 tmp737=tmp33*tmp623;
REAL8 tmp738=tmp736+tmp737;
REAL8 tmp739=0.5*tmp110*tmp113*tmp117*tmp33*tmp39*tmp738;
REAL8 tmp740=tmp726+tmp739;
REAL8 tmp717=2.*tmp500;
REAL8 tmp718=1.+tmp717;
REAL8 tmp214=(1.0/(tmp82*tmp82));
REAL8 tmp719=-(tmp110*tmp17*tmp222*tmp32*tmp33*tmp417*tmp420*tmp450*tmp47*tmp486*tmp50*tmp716*tmp718);
REAL8 tmp720=tmp186*tmp255*tmp259*tmp265;
REAL8 tmp721=1./sqrt(tmp720);
REAL8 tmp722=-2.*tmp11*tmp43*tmp82;
REAL8 tmp723=tmp513*tmp623;
REAL8 tmp724=tmp722+tmp723;
REAL8 tmp725=-0.5*tmp222*tmp32*tmp415*tmp420*tmp501*tmp522*tmp721*tmp724;
REAL8 tmp729=tmp169*tmp222*tmp417*tmp517*tmp728;
REAL8 tmp730=-(tmp116*tmp130*tmp17*tmp47*tmp50);
REAL8 tmp731=tmp169*tmp222*tmp728;
REAL8 tmp732=-(tmp116*tmp17*tmp47);
REAL8 tmp733=tmp110*tmp116*tmp17*tmp33*tmp47*tmp716;
REAL8 tmp734=tmp732+tmp733;
REAL8 tmp735=tmp130*tmp50*tmp734;
REAL8 tmp741=-(tmp169*tmp222*tmp740);
REAL8 tmp742=tmp731+tmp735+tmp741;
REAL8 tmp743=tmp500*tmp742;
REAL8 tmp744=tmp730+tmp743;
REAL8 tmp745=tmp417*tmp517*tmp744;
REAL8 tmp746=tmp222*tmp32*tmp415*tmp420*tmp522*tmp718*tmp740;
REAL8 tmp747=tmp729+tmp745+tmp746;
REAL8 tmp748=tmp417*tmp747;
REAL8 tmp749=tmp725+tmp748;
REAL8 tmp750=tmp513*tmp749;
REAL8 tmp751=tmp719+tmp750;
REAL8 tmp753=1/tmp506;
REAL8 tmp217=(1.0/(tmp43*tmp43));
REAL8 tmp676=2.*eta*tmp222*tmp43*tmp64*tmp72;
REAL8 tmp677=-2.*eta*tmp621*tmp64*tmp72;
REAL8 tmp678=-5.*coeffs->k5l*tmp59*x->data[2];
REAL8 tmp679=-5.*tmp59*tmp602*x->data[2];
REAL8 tmp680=tmp678+tmp679;
REAL8 tmp681=1.*tmp19*tmp680;
REAL8 tmp682=-(tmp59*tmp610*x->data[2]);
REAL8 tmp683=tmp681+tmp682;
REAL8 tmp684=1.*tmp19*tmp683;
REAL8 tmp685=-(tmp59*tmp612*x->data[2]);
REAL8 tmp686=tmp684+tmp685;
REAL8 tmp687=1.*tmp19*tmp686;
REAL8 tmp688=-(tmp59*tmp614*x->data[2]);
REAL8 tmp689=tmp687+tmp688;
REAL8 tmp690=-(eta*tmp43*tmp605*tmp689);
REAL8 tmp691=-5.*tmp51*tmp602*x->data[2];
REAL8 tmp692=tmp53+tmp54+tmp56+tmp58+tmp60+tmp691;
REAL8 tmp693=(1.0/(tmp604*tmp604));
REAL8 tmp694=eta*tmp43*tmp616*tmp692*tmp693;
REAL8 tmp695=-(eta*tmp605*tmp616*tmp86);
REAL8 tmp696=2.*tmp17*tmp59*tmp82*x->data[2];
REAL8 tmp697=2.*tmp222*tmp82*tmp86;
REAL8 tmp698=2.*tmp19*tmp43*tmp82*x->data[2];
REAL8 tmp699=tmp676+tmp677+tmp690+tmp694+tmp695+tmp696+tmp697+tmp698;
REAL8 tmp763=-(eta*tmp11*tmp43*tmp64*tmp72);
REAL8 tmp764=-(tmp11*tmp82*tmp86);
REAL8 tmp765=-2.*tmp43*tmp82*x->data[2];
REAL8 tmp766=tmp171+tmp763+tmp764+tmp765;
REAL8 tmp803=-(tmp173*tmp174*tmp222);
REAL8 tmp804=tmp116*tmp19*x->data[2];
REAL8 tmp805=(1.0/sqrt(tmp512*tmp512*tmp512));
REAL8 tmp806=0.5*tmp570*tmp805;
REAL8 tmp807=tmp803+tmp804+tmp806;
REAL8 tmp827=-4.*eta*tmp224*tmp43*tmp64*tmp72;
REAL8 tmp828=-4.*tmp224*tmp82*tmp86;
REAL8 tmp829=-12.*tmp222*tmp43*tmp82*x->data[2];
REAL8 tmp830=tmp33*tmp699;
REAL8 tmp831=2.*tmp623*x->data[2];
REAL8 tmp832=tmp827+tmp828+tmp829+tmp830+tmp831;
REAL8 tmp833=0.5*tmp110*tmp113*tmp117*tmp33*tmp39*tmp832;
REAL8 tmp834=-0.5*tmp113*tmp117*tmp33*tmp39*tmp738*tmp89*tmp93;
REAL8 tmp835=-0.5*eta*tmp110*tmp113*tmp214*tmp33*tmp39*tmp64*tmp72*tmp738;
REAL8 tmp836=-0.5*tmp110*tmp117*tmp217*tmp33*tmp39*tmp738*tmp86;
REAL8 tmp837=1.*tmp110*tmp113*tmp117*tmp39*tmp738*x->data[2];
REAL8 tmp838=-(tmp110*tmp113*tmp117*tmp33*tmp57*tmp738*x->data[2]);
REAL8 tmp839=tmp803+tmp804+tmp833+tmp834+tmp835+tmp836+tmp837+tmp838;
REAL8 tmp1=s1Vec->data[0]*s1Vec->data[0];
REAL8 tmp2=s1Vec->data[1]*s1Vec->data[1];
REAL8 tmp3=s1Vec->data[2]*s1Vec->data[2];
REAL8 tmp4=s2Vec->data[0]*s2Vec->data[0];
REAL8 tmp5=s2Vec->data[1]*s2Vec->data[1];
REAL8 tmp6=s2Vec->data[2]*s2Vec->data[2];
REAL8 tmp7=tmp1+tmp2+tmp3+tmp4+tmp5+tmp6;
REAL8 tmp119=1./sqrt(tmp118);
REAL8 tmp212=sqrt(tmp208);
REAL8 tmp490=tmp430*tmp430;
REAL8 tmp491=tmp436*tmp436;
REAL8 tmp492=tmp442*tmp442;
REAL8 tmp493=tmp486*tmp486;
REAL8 tmp494=-3.*tmp493;
REAL8 tmp495=tmp490+tmp491+tmp492+tmp494;
REAL8 d001000=(1.*eta*(2.*tmp109*tmp11*tmp110*tmp18+2.*tmp110*tmp18*tmp222*tmp413+tmp110*tmp135*tmp222*tmp32*tmp413*tmp415*tmp422*tmp424+tmp109*tmp110*tmp135*tmp222*tmp415*tmp422*tmp424*tmp444-0.5*tmp110*tmp135*tmp222*tmp237*tmp32*tmp415*tmp422*tmp444*tmp447+0.5*tmp110*tmp135*tmp173*tmp222*tmp32*tmp422*tmp424*tmp444*tmp450+tmp110*tmp135*tmp222*tmp32*tmp415*tmp424*tmp444*(0.5*tmp173*tmp420*tmp450-0.5*tmp455*tmp459+0.5*tmp415*tmp462*tmp468)-0.5*(2.*tmp345*tmp430+2.*tmp378*tmp436+2.*tmp411*tmp442-6.*tmp482*tmp486)*tmp59-tmp113*tmp117*tmp135*tmp173*tmp174*tmp39*tmp420*tmp751*tmp753+0.5*tmp113*tmp116*tmp117*tmp135*tmp39*tmp462*tmp468*tmp751*tmp753-eta*tmp113*tmp116*tmp135*tmp214*tmp39*tmp420*tmp64*tmp72*tmp751*tmp753+2.*tmp110*tmp177*tmp222*tmp32*tmp415*tmp422*tmp424*tmp444*tmp47*tmp78+2.*tmp113*tmp116*tmp117*tmp177*tmp39*tmp420*tmp47*tmp751*tmp753*tmp78-tmp116*tmp117*tmp135*tmp217*tmp39*tmp420*tmp751*tmp753*tmp86-tmp135*tmp424*tmp455*tmp462*tmp498*tmp513*tmp528*tmp530*tmp627*tmp671*tmp89+4.*tmp224*tmp424*tmp43*tmp455*tmp462*tmp47*tmp497*tmp498*tmp526*tmp528*tmp530*tmp82*tmp89-0.25*tmp135*tmp237*tmp455*tmp462*tmp498*tmp499*tmp502*tmp513*tmp627*tmp671*tmp93-0.25*tmp135*tmp237*tmp447*tmp455*tmp462*tmp498*tmp513*tmp528*tmp627*tmp671*tmp93-0.75*tmp135*tmp173*tmp424*tmp455*tmp462*tmp513*tmp528*tmp533*tmp627*tmp671*tmp93-0.25*tmp135*tmp424*tmp459*tmp462*tmp498*tmp513*tmp528*tmp538*tmp627*tmp671*tmp93-0.25*tmp135*tmp424*tmp455*tmp468*tmp498*tmp513*tmp528*tmp540*tmp627*tmp671*tmp93+0.25*tmp135*tmp424*tmp455*tmp462*tmp498*tmp528*tmp565*tmp570*tmp627*tmp671*tmp93-2.*eta*tmp224*tmp424*tmp43*tmp455*tmp462*tmp47*tmp497*tmp498*tmp526*tmp528*tmp64*tmp72*tmp93+1.*tmp177*tmp424*tmp455*tmp462*tmp47*tmp498*tmp513*tmp528*tmp627*tmp671*tmp78*tmp93+1.*tmp224*tmp237*tmp43*tmp455*tmp462*tmp47*tmp497*tmp498*tmp499*tmp502*tmp526*tmp82*tmp93+1.*tmp224*tmp237*tmp43*tmp447*tmp455*tmp462*tmp47*tmp497*tmp498*tmp526*tmp528*tmp82*tmp93+3.*tmp173*tmp224*tmp424*tmp43*tmp455*tmp462*tmp47*tmp497*tmp526*tmp528*tmp533*tmp82*tmp93+1.*tmp224*tmp424*tmp43*tmp459*tmp462*tmp47*tmp497*tmp498*tmp526*tmp528*tmp538*tmp82*tmp93+1.*tmp224*tmp424*tmp43*tmp455*tmp468*tmp47*tmp497*tmp498*tmp526*tmp528*tmp540*tmp82*tmp93-2.*tmp224*tmp424*tmp43*tmp455*tmp462*tmp497*tmp498*tmp526*tmp528*tmp78*tmp82*tmp93-2.*tmp224*tmp424*tmp455*tmp462*tmp47*tmp497*tmp498*tmp526*tmp528*tmp82*tmp86*tmp93-2.*tmp11*tmp18*tmp32*tmp89*tmp93-2.*tmp18*tmp222*tmp444*tmp89*tmp93-tmp135*tmp222*tmp32*tmp415*tmp422*tmp424*tmp444*tmp89*tmp93+4.*tmp110*tmp18*tmp32*x->data[2]+2.*tmp110*tmp18*tmp19*tmp444*x->data[2]+tmp110*tmp135*tmp19*tmp32*tmp415*tmp422*tmp424*tmp444*x->data[2]+1.5*tmp495*tmp55*x->data[2]-4.*coeffs->dheffSSv2*eta*tmp12*tmp7*x->data[2]-2.*tmp113*tmp116*tmp117*tmp135*tmp420*tmp57*tmp751*tmp753*x->data[2]-6.*tmp222*tmp424*tmp43*tmp455*tmp462*tmp47*tmp497*tmp498*tmp526*tmp528*tmp82*tmp93*x->data[2]+0.5*tmp135*tmp424*tmp455*tmp462*tmp498*tmp513*tmp528*tmp671*tmp93*(2.*tmp18*tmp89-2.*tmp18*tmp19*tmp625*x->data[2]-2.*tmp18*tmp222*(-(tmp17*tmp50*tmp699)+2.*tmp17*tmp47*tmp623*tmp78+8.*tmp222*x->data[2]+4.*tmp19*tmp33*x->data[2]))+0.5*tmp135*tmp424*tmp455*tmp462*tmp498*tmp513*tmp528*tmp627*tmp93*(-(tmp109*tmp11*tmp169*tmp415*tmp417*tmp420*tmp517)-tmp11*tmp156*tmp32*tmp415*tmp417*tmp420*tmp517-0.5*tmp11*tmp169*tmp173*tmp32*tmp417*tmp420*tmp450*tmp517-0.5*tmp11*tmp169*tmp32*tmp415*tmp420*tmp455*tmp459*tmp517-0.5*tmp11*tmp169*tmp32*tmp415*tmp417*tmp462*tmp468*tmp517-tmp11*tmp169*tmp32*tmp415*tmp417*tmp420*tmp578+eta*tmp110*tmp197*tmp265*tmp43*tmp503*tmp522*tmp64*tmp72+eta*tmp11*tmp43*tmp50*tmp64*tmp661*tmp72+2.*tmp110*tmp115*tmp173*tmp197*tmp265*tmp43*tmp522*tmp82+2.*tmp109*tmp110*tmp265*tmp32*tmp43*tmp503*tmp522*tmp82+tmp110*tmp197*tmp265*tmp43*tmp503*tmp590*tmp82-2.*tmp11*tmp43*tmp47*tmp661*tmp78*tmp82+tmp110*tmp197*tmp265*tmp503*tmp522*tmp82*tmp86+tmp11*tmp50*tmp661*tmp82*tmp86-tmp197*tmp265*tmp43*tmp503*tmp522*tmp82*tmp89*tmp93-2.*tmp169*tmp32*tmp415*tmp417*tmp420*tmp517*x->data[2]+4.*tmp11*tmp110*tmp197*tmp43*tmp503*tmp522*tmp82*x->data[2]+2.*tmp43*tmp50*tmp661*tmp82*x->data[2]+tmp11*tmp43*tmp50*tmp82*(tmp130*tmp169*tmp222*tmp482*tmp513+tmp130*tmp156*tmp222*tmp486*tmp513+tmp126*tmp169*tmp222*tmp486*tmp513+tmp173*tmp506*tmp522+tmp115*tmp522*tmp552+0.5*tmp130*tmp169*tmp222*tmp486*tmp565*tmp570+tmp115*tmp506*tmp590-eta*tmp11*tmp180*tmp186*tmp43*tmp522*tmp64*tmp72-2.*tmp11*tmp126*tmp130*tmp186*tmp43*tmp522*tmp82-tmp11*tmp180*tmp190*tmp191*tmp43*tmp522*tmp82-tmp11*tmp180*tmp186*tmp43*tmp590*tmp82-tmp11*tmp180*tmp186*tmp522*tmp82*tmp86+tmp130*tmp169*tmp19*tmp486*tmp513*x->data[2]-2.*tmp180*tmp186*tmp43*tmp522*tmp82*x->data[2]))-2.*tmp224*tmp424*tmp43*tmp455*tmp462*tmp47*tmp497*tmp498*tmp528*tmp82*tmp93*(tmp482*tmp510-tmp126*tmp417*tmp513*tmp524-0.5*tmp130*tmp455*tmp459*tmp513*tmp524-0.5*tmp130*tmp417*tmp524*tmp565*tmp570-tmp130*tmp417*tmp513*(-(tmp109*tmp222*tmp415*tmp420*tmp517)-0.5*tmp173*tmp222*tmp32*tmp420*tmp450*tmp517-0.5*tmp222*tmp32*tmp415*tmp462*tmp468*tmp517+tmp156*tmp222*tmp417*tmp522+0.5*tmp169*tmp222*tmp455*tmp459*tmp522-tmp222*tmp32*tmp415*tmp420*tmp578+tmp169*tmp222*tmp417*tmp590-tmp19*tmp32*tmp415*tmp420*tmp517*x->data[2]+tmp169*tmp19*tmp417*tmp522*x->data[2])+tmp486*(-(eta*tmp110*tmp197*tmp265*tmp43*tmp503*tmp64*tmp72)+eta*tmp11*tmp43*tmp508*tmp64*tmp72-2.*tmp110*tmp115*tmp173*tmp197*tmp265*tmp43*tmp82-2.*tmp109*tmp110*tmp265*tmp32*tmp43*tmp503*tmp82-tmp110*tmp197*tmp265*tmp503*tmp82*tmp86+tmp11*tmp508*tmp82*tmp86+tmp197*tmp265*tmp43*tmp503*tmp82*tmp89*tmp93-4.*tmp11*tmp110*tmp197*tmp43*tmp503*tmp82*x->data[2]+2.*tmp43*tmp508*tmp82*x->data[2]+tmp11*tmp43*tmp82*(2.*tmp11*tmp156*tmp169-tmp173*tmp50*tmp506-tmp115*tmp50*tmp552+2.*tmp115*tmp47*tmp506*tmp78+2.*tmp175*x->data[2])))-0.5*tmp212*(tmp113*tmp116*tmp117*tmp39*tmp89-tmp113*tmp117*tmp173*tmp174*tmp39*tmp92-eta*tmp113*tmp116*tmp214*tmp39*tmp64*tmp72*tmp92-tmp116*tmp117*tmp217*tmp39*tmp86*tmp92-2.*tmp113*tmp116*tmp117*tmp57*tmp92*x->data[2])*(1.0/sqrt(tmp118*tmp118*tmp118))-tmp113*tmp116*tmp117*tmp135*tmp39*tmp420*tmp552*tmp751*(1.0/(tmp506*tmp506))+tmp113*tmp116*tmp117*tmp135*tmp39*tmp420*tmp753*(-(tmp110*tmp17*tmp222*tmp237*tmp32*tmp33*tmp417*tmp420*tmp424*tmp450*tmp47*tmp486*tmp50*tmp716)-tmp110*tmp17*tmp222*tmp32*tmp33*tmp417*tmp420*tmp450*tmp47*tmp482*tmp50*tmp716*tmp718-tmp109*tmp110*tmp17*tmp222*tmp33*tmp417*tmp420*tmp450*tmp47*tmp486*tmp50*tmp716*tmp718-0.5*tmp110*tmp17*tmp222*tmp32*tmp33*tmp420*tmp450*tmp455*tmp459*tmp47*tmp486*tmp50*tmp716*tmp718-0.5*tmp110*tmp17*tmp222*tmp32*tmp33*tmp417*tmp450*tmp462*tmp468*tmp47*tmp486*tmp50*tmp716*tmp718+0.5*tmp110*tmp17*tmp173*tmp222*tmp32*tmp33*tmp417*tmp420*tmp47*tmp486*tmp498*tmp50*tmp716*tmp718+0.5*tmp565*tmp570*tmp749-tmp110*tmp17*tmp222*tmp32*tmp33*tmp417*tmp420*tmp450*tmp47*tmp486*tmp50*tmp718*tmp766+2.*tmp110*tmp17*tmp222*tmp32*tmp33*tmp417*tmp420*tmp450*tmp48*tmp486*tmp716*tmp718*tmp78-tmp110*tmp17*tmp222*tmp32*tmp33*tmp417*tmp420*tmp450*tmp486*tmp50*tmp716*tmp718*tmp78+tmp17*tmp222*tmp32*tmp33*tmp417*tmp420*tmp450*tmp47*tmp486*tmp50*tmp716*tmp718*tmp89*tmp93-2.*tmp110*tmp17*tmp222*tmp32*tmp417*tmp420*tmp450*tmp47*tmp486*tmp50*tmp716*tmp718*x->data[2]-tmp110*tmp17*tmp19*tmp32*tmp33*tmp417*tmp420*tmp450*tmp47*tmp486*tmp50*tmp716*tmp718*x->data[2]+tmp513*(-0.25*tmp222*tmp237*tmp32*tmp415*tmp420*tmp424*tmp522*tmp721*tmp724-0.5*tmp109*tmp222*tmp415*tmp420*tmp501*tmp522*tmp721*tmp724-0.25*tmp173*tmp222*tmp32*tmp420*tmp450*tmp501*tmp522*tmp721*tmp724-0.25*tmp222*tmp32*tmp415*tmp462*tmp468*tmp501*tmp522*tmp721*tmp724-0.5*tmp222*tmp32*tmp415*tmp420*tmp501*tmp590*tmp721*tmp724+0.5*tmp455*tmp459*tmp747-0.5*tmp19*tmp32*tmp415*tmp420*tmp501*tmp522*tmp721*tmp724*x->data[2]-0.5*tmp222*tmp32*tmp415*tmp420*tmp501*tmp522*tmp721*(0.5*tmp565*tmp570*tmp623+tmp513*tmp699-2.*eta*tmp11*tmp43*tmp64*tmp72-2.*tmp11*tmp82*tmp86-4.*tmp43*tmp82*x->data[2])+tmp417*(tmp156*tmp222*tmp417*tmp517*tmp728+0.5*tmp169*tmp222*tmp455*tmp459*tmp517*tmp728+tmp169*tmp222*tmp417*tmp578*tmp728+1.*tmp222*tmp237*tmp32*tmp415*tmp420*tmp424*tmp522*tmp740+tmp109*tmp222*tmp415*tmp420*tmp522*tmp718*tmp740+0.5*tmp173*tmp222*tmp32*tmp420*tmp450*tmp522*tmp718*tmp740+0.5*tmp222*tmp32*tmp415*tmp462*tmp468*tmp522*tmp718*tmp740+tmp222*tmp32*tmp415*tmp420*tmp590*tmp718*tmp740+0.5*tmp455*tmp459*tmp517*tmp744+tmp417*tmp578*tmp744+tmp169*tmp222*tmp417*tmp517*tmp807+tmp222*tmp32*tmp415*tmp420*tmp522*tmp718*tmp839+tmp169*tmp19*tmp417*tmp517*tmp728*x->data[2]+tmp19*tmp32*tmp415*tmp420*tmp522*tmp718*tmp740*x->data[2]+tmp417*tmp517*(-(tmp116*tmp126*tmp17*tmp47*tmp50)+tmp130*tmp17*tmp173*tmp174*tmp47*tmp50+0.5*tmp237*tmp424*tmp742+2.*tmp116*tmp130*tmp17*tmp48*tmp78-tmp116*tmp130*tmp17*tmp50*tmp78+tmp500*(tmp156*tmp222*tmp728+tmp126*tmp50*tmp734-tmp156*tmp222*tmp740-2.*tmp130*tmp47*tmp734*tmp78+tmp169*tmp222*tmp807-tmp169*tmp222*tmp839+tmp169*tmp19*tmp728*x->data[2]-tmp169*tmp19*tmp740*x->data[2]+tmp130*tmp50*(tmp17*tmp173*tmp174*tmp47-tmp110*tmp17*tmp173*tmp174*tmp33*tmp47*tmp716+tmp110*tmp116*tmp17*tmp33*tmp47*tmp766-tmp116*tmp17*tmp78+tmp110*tmp116*tmp17*tmp33*tmp716*tmp78-tmp116*tmp17*tmp33*tmp47*tmp716*tmp89*tmp93+2.*tmp110*tmp116*tmp17*tmp47*tmp716*x->data[2]))))+0.25*tmp222*tmp32*tmp415*tmp420*tmp501*tmp522*tmp724*(tmp190*tmp191*tmp255*tmp259*tmp265+2.*eta*tmp186*tmp255*tmp265*tmp64*tmp72*tmp82+2.*tmp186*tmp259*tmp265*tmp43*tmp86+4.*tmp11*tmp186*tmp255*tmp259*x->data[2])*(1.0/sqrt(tmp720*tmp720*tmp720))))+(0.5*tmp119*(tmp170+tmp176+tmp178+tmp179+tmp187+tmp192+tmp193+tmp194+tmp195+tmp196+tmp198+tmp199+tmp200+tmp201+tmp202+8.*eta*tmp121*tmp126*tmp131*tmp39-4.*eta*tmp121*tmp133*tmp57*x->data[2]))/sqrt(tmp208)))/sqrt(1.+2.*eta*(-1.+tmp119*tmp212+2.*tmp11*tmp110*tmp18*tmp32+2.*tmp110*tmp18*tmp222*tmp444+tmp110*tmp135*tmp222*tmp32*tmp415*tmp422*tmp424*tmp444-0.5*tmp495*tmp59+coeffs->dheffSSv2*eta*tmp57*tmp7+tmp113*tmp116*tmp117*tmp135*tmp39*tmp420*tmp751*tmp753+0.5*tmp135*tmp424*tmp455*tmp462*tmp498*tmp513*tmp528*tmp627*tmp671*tmp93-2.*tmp224*tmp424*tmp43*tmp455*tmp462*tmp47*tmp497*tmp498*tmp526*tmp528*tmp82*tmp93));
