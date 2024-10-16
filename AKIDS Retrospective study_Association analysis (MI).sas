/*****************************/
/* Table2 Cox Association analysis */
/*****************************/

/**** 1. Cause-specific hazards models ****/
%macro Cox(AKIgroup, Day, Event, Output);
/*Model 3*/
proc phreg data=lib.Population_79838_MI_model3;
by _imp;
class &AKIgroup. (ref='1') gender(ref='0') MI_Smoking(ref='0') dm_define(ref='0') htn_define(ref='0') Chf(ref='0') Cvd(ref='0') 
        Dementia(ref='0')  Renal(ref='0') NSAID_b90d(ref='0') Contrast_b90d(ref='0') ACEI_ARBs_drug_b90d(ref='0')
		Antimicrobial_b90d(ref='0') Diuretics_b90d(ref='0') ; 
model &day.*&event.(0)=&AKIgroup. age gender MI_BMI MI_Smoking dm_define htn_define Chf Cvd Dementia Renal 
										NSAID_b90d Contrast_b90d ACEI_ARBs_drug_b90d Antimicrobial_b90d Diuretics_b90d 
										 MI_Hb_b1y min_SCr /  risklimits;
ods output ParameterEstimates=Model3;
run;

proc sort data=Model3;by parameter ClassVal0 _imp;run;
proc mianalyze data=Model3;
	by parameter ClassVal0;
	modeleffects estimate;
	stderr stderr;
	ods output ParameterEstimates = Model3a;
run;

data Cox_Model3_&output.;  set Model3a;
	Log_HR_comb=Estimate;
	HazardRatio=exp(Estimate);
	HRLowerCL=exp(LCLMean);
	HRUpperCL=exp(UCLMean);
	keep Parameter ClassVal0 Log_HR_comb HazardRatio HRLowerCL HRUpperCL Probt;
	rename Probt=HR_pval_comb;
run;

data Cox_Model3_&output._HR; set Cox_Model3_&output.;
HR=left(trim(round(HazardRatio,0.1))||" ("||trim(left(round(HRLowerCL,0.1)))||"-"||trim(left(round(HRUpperCL,0.1)))||")"||sign);
run;

%mend;
/*AKI_type_merge(two groups): Non-AKDOPT and AKDOPT*/
%Cox(AKI_type_merge, Lossfollow_CKO_FT, CKO_a1y, CKO_a1y_merge); 
%Cox(AKI_type_merge, Lossfollow_death_FT, death_a1y, death_a1y_merge); 
%Cox(AKI_type_merge, Lossfollow_De_Novo_FT, De_Novo_CKD_2, De_Novo_CKD_2_merge); 

/*AKI_type(three groups): Non-AKDOPT and Stable AKDOPT and Deteriorating AKDOPT*/
%Cox(AKI_type, Lossfollow_CKO_FT, CKO_a1y, CKO_a1y); 
%Cox(AKI_type, Lossfollow_death_FT, death_a1y, death_a1y); 
%Cox(AKI_type, Lossfollow_De_Novo_FT, De_Novo_CKD_2, De_Novo_CKD_2); 


/**** 2. Estimating Absolute Risk ****/
%macro MI;
%do i=1 %to 20;
	data MI_&i.;set lib.Population_79838_mi_model3; 
	if _imp=&i. then output MI_&i.;
	keep PatientNo AKI_type_merge AKI_type
				Lossfollow_CKO_FT CKO_a1y  
				Lossfollow_death_FT death_a1y
				Lossfollow_De_Novo_FT De_Novo_CKD_2
				age gender MI_BMI MI_Smoking dm_define htn_define Chf Cvd Dementia Renal
				NSAID_b90d Contrast_b90d ACEI_ARBs_drug_b90d Antimicrobial_b90d Diuretics_b90d
				MI_Hb_b1y min_SCr ;
	run; 
%end;
%mend; 
%MI;
 
%macro pop ;
%do i=1 %to 20  ;
		/*Assume everyone is Non-AKI*/
		data AA;set MI_&i.(drop=AKI_type_merge AKI_type);AKI_type_merge=1;run; 
		/*Assume everyone is AKI*/
		data BB;set MI_&i.(drop=AKI_type_merge AKI_type);AKI_type_merge=2;run; 
		/*Merge (AKI_type_merge)*/
		data MI_&i._Population1;set AA BB; run; 
		/*Clean*/
		proc delete data=AA BB;run;

		/*Assume everyone is Non-AKI*/
		data AA;set MI_&i.(drop=AKI_type AKI_type_merge);AKI_type=1;run; 
		/*Assume everyone is S-AKI*/
		data BB;set MI_&i.(drop=AKI_type AKI_type_merge);AKI_type=2;run; 
		/*Assume everyone is D-AKI*/
		data CC;set MI_&i.(drop=AKI_type AKI_type_merge);AKI_type=3;run; 
		/*Merge (AKI_type)*/
		data MI_&i._Population2;set AA BB CC; run; 
		/*Clean*/
		proc delete data=AA BB CC;run;
%end;
%mend;
%pop ; 


/*    Fit Cox proportional hazards regression model.   */
  %macro model ( follow, outcome );
 %do i=1 %to 20  ;
/* AKI type merge (Non-AKI and AKI) */
proc phreg data= MI_&i.;
	class AKI_type_merge(ref='1') gender(ref='0') MI_Smoking(ref='0') dm_define(ref='0') htn_define(ref='0') Chf(ref='0') Cvd(ref='0') 
	        Dementia(ref='0')  Renal(ref='0') NSAID_b90d(ref='0') Contrast_b90d(ref='0') ACEI_ARBs_drug_b90d(ref='0')
			Antimicrobial_b90d(ref='0') Diuretics_b90d(ref='0')  ; 
	model &follow.*&outcome.(0)=AKI_type_merge age gender MI_BMI MI_Smoking dm_define htn_define Chf Cvd Dementia Renal 
											NSAID_b90d Contrast_b90d ACEI_ARBs_drug_b90d Antimicrobial_b90d Diuretics_b90d 
											MI_Hb_b1y min_SCr ;
		baseline /*Establishes a baseline risk function*/ 
		out=&outcome._2group  covariates=MI_&i._Population1 
		survival=survival/nomean /*Nomean does not need to calculate the average predicted survival curve*/; 
run;

proc sort data=&outcome._2group;by PatientNo AKI_type_merge &follow.;run; 
data &outcome._2group;set &outcome._2group;
	by PatientNo AKI_type_merge Lossfollow_CKO_FT;
	if last.AKI_type_merge;
	event_risk=1 - survival; *Calculate risk probability; 
run;

proc sort data=&outcome._2group; by AKI_type_merge; run;
proc means mean data=&outcome._2group noprint;
	var event_risk;
	by AKI_type_merge;
	output out=pop_risk1 mean=pop_risk1; /* Compute the mean probability of outcome  */
run;
proc transpose data=pop_risk1 out=pop_risk1;var pop_risk1;run;
data temp.MI_&i._&outcome._2group;
set pop_risk1;
	AKI_type_1=col1;	AKI_type_2=col2;
	Risk_difference=AKI_type_2 - AKI_type_1;	/* Risk difference */
	Outcome="&outcome.";
	MI_Group=&i.;
keep AKI_type_1 AKI_type_2 Risk_difference Outcome MI_Group;
run;

/*Clean*/
proc delete data=&outcome._2group pop_risk1;run;
/********************************************************************/
/* AKI type (Non-AKI and S-AKI and D-AKI)*/
proc phreg data= MI_&i.;
	class AKI_type(ref='1') gender(ref='0') MI_Smoking(ref='0') dm_define(ref='0') htn_define(ref='0') Chf(ref='0') Cvd(ref='0') 
	        Dementia(ref='0')  Renal(ref='0') NSAID_b90d(ref='0') Contrast_b90d(ref='0') ACEI_ARBs_drug_b90d(ref='0')
			Antimicrobial_b90d(ref='0') Diuretics_b90d(ref='0')  ; 
	model &follow.*&outcome.(0)=AKI_type age gender MI_BMI MI_Smoking dm_define htn_define Chf Cvd Dementia Renal 
											NSAID_b90d Contrast_b90d ACEI_ARBs_drug_b90d Antimicrobial_b90d Diuretics_b90d 
											MI_Hb_b1y min_SCr ;
		baseline /*Establishes a baseline risk function*/ 
		out=&outcome._3group  covariates=MI_&i._Population2 
		survival=survival/nomean /*Nomean does not need to calculate the average predicted survival curve*/; 
run;

proc sort data=&outcome._3group;by PatientNo AKI_type &follow.;run; 
data &outcome._3group;set &outcome._3group;
	by PatientNo AKI_type Lossfollow_CKO_FT;
	if last.AKI_type;
	event_risk=1 - survival; *Calculate risk probability; 
run;
proc sort data=&outcome._3group; by AKI_type; run;
proc means mean data=&outcome._3group noprint;
	var event_risk;
	by AKI_type;
	output out=pop_risk2 mean=pop_risk2; /* Compute the mean probability of outcome  */
run;
proc transpose data=pop_risk2 out=pop_risk2; var pop_risk2; run;
data  temp.MI_&i._&outcome._3group;
set pop_risk2;
	AKI_type_1=col1;
	AKI_type_2=col2;
	AKI_type_3=col3;
/*Non-AKI vs S-AKI*/
	Risk_difference=AKI_type_2 - AKI_type_1;	/* Risk difference */
/*Non-AKI vs D-AKI*/
	Risk_difference_2=AKI_type_3 - AKI_type_1;	/* Risk difference */
	Outcome="&outcome.";
	MI_Group=&i.;
keep AKI_type_1 AKI_type_2 AKI_type_3 Risk_difference Risk_difference_2 Outcome MI_Group;
run;

/*Clean */
proc delete data=&outcome._3group pop_risk2;run; 
%end;
%mend;
%model( Lossfollow_CKO_FT, CKO_a1y ); 
%model(Lossfollow_death_FT, death_a1y);
%model(Lossfollow_De_Novo_FT, De_Novo_CKD_2);

/*Integrate the dataset*/
Data Absolute_risk_AKI_type_merge; set _null_ ;run; 
Data Absolute_risk_AKI_type; set _null_ ;length Outcome $30.;run; 
%macro Absolute  ;
 %do i=1 %to 20;
data Absolute_risk_AKI_type_merge; length Outcome $30.;
	Set Absolute_risk_AKI_type_merge 
	Temp.MI_&i._Cko_a1y_2group
	 Temp.MI_&i._Death_a1y_2group
	 Temp.MI_&i._de_novo_ckd_2_2group;
run;

data Absolute_risk_AKI_type; length Outcome $30.;
	Set Absolute_risk_AKI_type 
	Temp.MI_&i._Cko_a1y_3group 
	 Temp.MI_&i._Death_a1y_3group
	 Temp.MI_&i._de_novo_ckd_2_3group;
run;
%end;
%mend;
%Absolute;
/*Calculate absolute risk*/
 proc means data=Absolute_risk_aki_type_merge mean ;var Risk_difference Risk_difference_2 ; class Outcome;run; 
 proc means data=Absolute_risk_AKI_type mean ;var Risk_difference ; class Outcome;run; 
