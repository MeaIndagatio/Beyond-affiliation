****** Data preparation: 2023-24 Religious Landscape Study *******



**# Load file

use "path/to/data/RLS Public Use File_1_original.dta", clear



**# Save log of all data preparation

log using "path/to/data/Log file of data preparation", replace



**# Data Preparation

	* Make lowercase and add number labels to all variables

	rename *, lower
	numlabel _all, add 



	* Check dataset for variables with refusals and set refusals to missing

	foreach x of varlist _all {
		capture local lbl : value label `x'
		if "`lbl'" != "" {
			quietly levelsof `x', local(vals)
			foreach v of local vals {
				local t : label `lbl' `v'
				if strpos(lower("`t'"), "refus") | strpos(lower("`t'"), "know") | strpos(lower("`t'"), "not sure") | strpos(lower("`t'"), "missing") {
					di "`x'  `v'  `t'"
				}
			}
		}
	}


	foreach x of varlist happy satis_a satis_b reltrad frmrel ///
    relper spirper attndperrls attndonrls relimp memb god god2 hvn hll soul bibimp ///
    pray grace prac_a prac_b prac_c prac_d ///
    exp_a exp_b exp_c exp_d exp_e exp_f exp_g ///
    spiract_a spiract_b spiract_c spiract_d spiract_e spiract_f rtrt ///
    scrlcon1 spirworld2 secbel1 secbel2 lifedir lifespir ///
    chimprel_a chimprel_b chattend chedu1 chedu2 ///
    marital sprelimp party partyln ideo inc_sdt1 gender ///
    birthdecade educrec usgen {
    recode `x' 99=. 900000=.
}



	* Preparing the variables for analysis

	summarize happy satis_a satis_b reltrad frmrel relper spirper attndperrls attndonrls ///
	relimp memb god god2 hvn hll soul bibimp pray grace prac_a prac_b prac_c prac_d ///
	exp_a exp_b exp_c exp_d exp_e exp_f exp_g spiract_a spiract_b spiract_c spiract_d ///
	spiract_e spiract_f rtrt scrlcon1 spirworld2 secbel1 secbel2 lifedir lifespir ///
	chimprel_a chimprel_b chattend chedu1 chedu2 marital sprelimp party partyln ideo ///
	inc_sdt1 gender birthdecade educrec usgen
	


	/* 1. Identifier and survey design 
 	p_suid: respondent ID.
	weight: full-sample weight.
	repwt_001 to repwt_500: bootstrap replicate weights. */

	ren p_suid respid

	
* 2. LPA indicator items

	* Religious behaviour items
	
	codebook attndperrls attndonrls memb pray grace prac_a prac_b prac_c prac_d

	
	* Religious behaviour index

	// i. Reverse-code so higher = more

	recode attndperrls (1=6) (2=5) (3=4) (4=3) (5=2) (6=1), gen(attndperrlsr)
	la var attndperrlsr "How often do you attend religious services in person (rev)"

	recode attndonrls (1=6) (2=5) (3=4) (4=3) (5=2) (6=1), gen(attndonrlsr)
	la var attndonrlsr "How often do you watch or participate in religious services online (rev)"

	la def attend 6 "More than once a week" 5 "Once a week" 4 "Once or twice a month" ///
		3 "A few times a year" 2 "Seldom" 1 "Never"
	la val attndperrlsr attndonrlsr attend

	recode memb (2=1 "no") (1=2 "yes"), gen(synmemr)
	la var synmemr "Are you personally a member of a house of worship (rev)"

	recode pray (1=7 "Several times a day") (2=6 "Once a day") (3=5 "A few times a week") ///
		(4=4 "Once a week") (5=3 "A few times a month") (6=2 "Seldom") ///
		(7=1 "Never"), gen(prayfreqr)
	la var prayfreqr "Outside of attending religious services, how often do you pray? (rev)"

	recode grace (1=5 "Always") (2=4 "Often") (3=3 "Sometimes") (4=2 "Seldom") ///
		(5=1 "Never"), gen(gracefreqr)
	la var gracefreqr "How often do you say grace or pray before meals? (rev)"

	foreach x of varlist prac_a prac_b prac_c prac_d {
		recode `x' (1=5) (2=4) (3=3) (4=2) (5=1), gen(`x'r)
	}
	la var prac_ar "How often do you participate in prayer groups, scripture study groups... (rev)"
	la var prac_br "How often do you read scripture outside of religious services? (rev)"
	la var prac_cr "How often do you read other holy texts, devotions (rev)"
	la var prac_dr "How often do you share your views about religion with people from other (rev)"
	la def prac 5 "At least once a week" 4 "Once or twice a month" ///
		3 "Several times a year" 2 "Seldom" 1 "Never"
	la val prac_ar prac_br prac_cr prac_dr prac

	// ii. Check recodes (each table should be a diagonal)

	tab attndperrls attndperrlsr, mi
	tab attndonrls attndonrlsr, mi
	tab memb synmemr, mi
	tab pray prayfreqr, mi
	tab grace gracefreqr, mi
	tab prac_a prac_ar, mi
	tab prac_b prac_br, mi
	tab prac_c prac_cr, mi
	tab prac_d prac_dr, mi

	// iii. Build the index (complete cases only)

	egen rb_nmiss = rowmiss(attndperrlsr attndonrlsr synmemr prayfreqr gracefreqr prac_ar prac_br prac_cr prac_dr)
	egen religbehav = rowtotal(attndperrlsr attndonrlsr synmemr prayfreqr gracefreqr prac_ar prac_br prac_cr prac_dr) if rb_nmiss == 0
	la var religbehav "Religious behaviours scale (sum of 9 items)"

	// iv. Diagnostics

	tab rb_nmiss
	alpha attndperrlsr attndonrlsr synmemr prayfreqr gracefreqr prac_ar prac_br prac_cr prac_dr, item
	summarize religbehav, detail

	
	
	/* Religious behavior was measured with nine items from the 2023–24 Pew Research Center 
	Religious Landscape Study:

	* frequency of in-person attendance at religious services (1 = never to 6 = more than 
	once a week);
	* frequency of watching or participating in religious services online or on television 
	(1–6);
	membership in a house of worship (1 = no, 2 = yes);
	* frequency of prayer outside religious services (1 = never to 7 = several times a day);
	* frequency of saying grace before meals (1 = never to 5 = always);
	* frequency of participation in prayer, scripture study, or religious education groups 
	(1 = never to 5 = at least once a week);
	* frequency of reading scripture outside religious services (1–5);
	* frequency of reading other holy texts, devotions, or inspirational literature (1–5);
	* frequency of sharing views about religion with people from other religious 
	backgrounds (1–5).

	"Don't know" and refused responses were coded as missing. Each item was reverse-coded 
	so that higher values represented more frequent practice. The index was calculated as 
	the sum of the nine items, giving a possible range of 9 to 46. It was computed only for 
	respondents with complete data on all nine items (n = 35,522; 96.2% of the sample of 
	36,908).

	The index demonstrated high internal consistency (Cronbach's α = .92). Corrected 
	item–total correlations ranged from .61 (sharing religious views) to .85 (reading 
	scripture). Alpha if an item was deleted ranged from .90 to .92, indicating that no 
	item substantially weakened the scale. Scores covered the full possible range (M = 
	22.82, SD = 10.76, Mdn = 21), with a modest positive skew (skewness = 0.42). At least 
	5% of respondents obtained the minimum score of 9, indicating a floor at the absence of 
	religious practice. Reliability and descriptive statistics are unweighted. */
	
	
	
	* Religious belief index (decided to leave out secbel2)

	codebook relimp relper god god2 hvn hll bibimp 
	
	
	// i. Reverse-code so higher = more

	recode relimp (1=4 "Very important") (2=3 "Somewhat important") ///
		(3=2 "Not too important") (4=1 "Not at all important"), gen(relimpr)
	la var relimpr "How important is religion in your life? (rev)"

	recode relper (1=4 "Very religious") (2=3 "Somewhat religious") ///
		(3=2 "Not too religious") (4=1 "Not at all religious"), gen(relperr)
	la var relperr "To what extent do you consider yourself a religious person? (rev)"

	recode hvn (1=2 "yes") (2=1 "no"), gen(heavenr)
	la var heavenr "Do you think there is a heaven (rev)"

	recode hll (1=2 "yes") (2=1 "no"), gen(hellr)
	la var hellr "Do you think there is a hell (rev)"
	
	/// Afterlife belief (heaven and hell combined)
	tab heavenr hellr, mi
	gen afterlife = heavenr + hellr - 1
	la def afterlife 1 "Neither heaven nor hell" 2 "Heaven or hell, not both" ///
		3 "Both heaven and hell"
	la val afterlife afterlife
	la var afterlife "Belief in heaven and hell (1-3)"
	tab afterlife, mi
	
	recode secbel2 (1=2 "Necessary to believe in God to be moral") ///
		(2=1 "Possible to be moral without believing in God"), gen(secbel2r)
	la var secbel2r "In order to be good and moral (God necessary or not) (rev)"

	// ii. Certainty of belief in God (corrected: 1 = does not believe, 2-5 = believer's certainty)

	gen gdcertainr = .
	replace gdcertainr = 1 if god == 2
	replace gdcertainr = 2 if god == 1 & god2 == 4
	replace gdcertainr = 3 if god == 1 & god2 == 3
	replace gdcertainr = 4 if god == 1 & god2 == 2
	replace gdcertainr = 5 if god == 1 & god2 == 1
	la def gdcertainr 1 "Does not believe" 2 "Believes, not at all certain" ///
		3 "Believes, not too certain" 4 "Believes, fairly certain" ///
		5 "Believes, absolutely certain"
	la val gdcertainr gdcertainr
	la var gdcertainr "Belief in God and certainty of belief (1-5)"

	// iii. Check recodes

	tab relimp relimpr, mi
	tab relper relperr, mi
	tab hvn heavenr, mi
	tab hll hellr, mi
	tab secbel2 secbel2r, mi
	tab god2 god, mi
	tab gdcertainr, mi

	// iv. Build the index (complete cases only)

	egen rl_nmiss = rowmiss(relimpr relperr gdcertainr afterlife)
	egen religbelif = rowtotal(relimpr relperr gdcertainr afterlife) if rl_nmiss == 0
	la var religbelif "Religious belief scale (sum of 4 items)"

	// v. Diagnostics

	tab rl_nmiss
	alpha relimpr relperr gdcertainr heavenr hellr, item
	summarize religbelif, detail
	


	* Spiritual behaviour index

	codebook spiract_a spiract_b spiract_c spiract_d spiract_e spiract_f rtrt

	

	// i. Spiritual behavior index (original items, complete cases only)

	foreach x of varlist spiract_a spiract_b spiract_c spiract_d spiract_e spiract_f {
		recode `x' (1=5) (2=4) (3=3) (4=2) (5=1), gen(`x'r)
	}
	la var spiract_ar "How often do you visit a nature spot for spiritual reasons? (rev)"
	la var spiract_br "How often do you listen to music for spiritual reasons? (rev)"
	la var spiract_cr "How often do you exercise for spiritual reasons? (rev)"
	la var spiract_dr "How often do you look inward or center yourself for spiritual reasons? (rev)"
	la var spiract_er "How often do you practice yoga for spiritual reasons? (rev)"
	la var spiract_fr "How often do you meditate for spiritual reasons? (rev)"
	la def spiract 5 "Do this daily for spiritual reasons" 4 "Weekly" 3 "Monthly" ///
		2 "Less often" 1 "Never do this for spiritual reasons"
	la val spiract_ar spiract_br spiract_cr spiract_dr spiract_er spiract_fr spiract

	recode rtrt (1=5 "More than once a year") (2=4 "About once a year") ///
		(3=3 "Once every few years") ///
		(4=2 "Less often") (5=1 "Never"), gen(retreatr)
	la var retreatr "How often do you go on spiritual retreats? (rev)"

	// ii. Check recodes (each table should be a diagonal)

	tab spiract_a spiract_ar, mi
	tab spiract_b spiract_br, mi
	tab spiract_c spiract_cr, mi
	tab spiract_d spiract_dr, mi
	tab spiract_e spiract_er, mi
	tab spiract_f spiract_fr, mi
	tab rtrt retreatr, mi

	// iii. Build the index (complete cases only) reanalysis 261021 dropped yoga

	egen sb_nmiss = rowmiss(spiract_ar spiract_br spiract_cr spiract_dr spiract_fr)
	egen spirbehav = rowtotal(spiract_ar spiract_br spiract_cr spiract_dr spiract_fr) if sb_nmiss == 0
	la var spirbehav "Spiritual behaviours scale (sum of 5 items)"

	// iv. Diagnostics

	tab sb_nmiss
	alpha spiract_ar spiract_br spiract_cr spiract_dr spiract_fr, item
	alpha spiract_ar spiract_br spiract_cr spiract_dr spiract_fr retreatr, item
	summarize spirbehav, detail
	
	
	
	* Spiritual belief index (5 items)

	codebook spirper soul spirworld2 secbel1 exp_a exp_b exp_c exp_d exp_e exp_f exp_g
	
	// i. Reverse-code so higher = more

	recode spirper (1=4 "Very spiritual") (2=3 "Somewhat spiritual") ///
		(3=2 "Not too spiritual") ///
		(4=1 "Not at all spiritual"), gen(spirperr)
	la var spirperr "To what extent do you consider yourself a spiritual person? (rev)"

	recode soul (1=2 "yes") (2=1 "no"), gen(soulr)
	la var soulr "Do you believe people have a soul or spirit (rev)"

	* spirworld2 is used as coded (1 = natural world is all there is, 2 = something spiritual beyond)

	recode exp_a (1=5 "At least once a week") (2=4 "Once or twice a month") ///
		(3=3 "Several times a year") ///
		(4=2 "Seldom") (5=1 "Never"), gen(exp_ar)
	la var exp_ar "How often do you feel a deep sense of spiritual peace and well-being? (rev)"

	recode exp_d (1=5 "At least once a week") (2=4 "Once or twice a month") ///
		(3=3 "Several times a year") ///
		(4=2 "Seldom") (5=1 "Never"), gen(exp_dr)
	la var exp_dr "How often do you feel the presence of something from beyond this world? (rev)"

	// ii. Check recodes

	tab spirper spirperr, mi
	tab soul soulr, mi
	tab spirworld2, mi
	tab exp_a exp_ar, mi
	tab exp_d exp_dr, mi

	// iii. Build the index (complete cases only)
		
	egen sl_nmiss = rowmiss(spirperr soulr spirworld2 exp_ar exp_dr)
	egen spirbelif = rowtotal(spirperr soulr spirworld2 exp_ar exp_dr) if sl_nmiss == 0
	la var spirbelif "Spiritual belief and experience scale (sum of 5 items)"

	// iv. Diagnostics

	tab sl_nmiss
	alpha spirperr soulr spirworld2 exp_ar exp_dr, item
	summarize spirbelif, detail

	
	/* Spiritual belief 
	Spiritual belief and experience. Spiritual belief and experience was measured with five 
	items:

	* the extent to which the respondent considers themselves a spiritual person (1 = not 
	at all spiritual to 4 = very spiritual);
	* belief that people have a soul or spirit in addition to their physical body (1 = no, 
	2 = yes);
	* belief that "there is something spiritual beyond the natural world, even if we cannot 
	see it" (2), as opposed to "the natural world is all there is" (1);
	* frequency of feeling a deep sense of spiritual peace and well-being (1 = never to 5 = 
	at least once a week);
	* frequency of feeling the presence of something from beyond this world (1 = never to 5 
	= at least once a week).

	"Don't know" and refused responses were coded as missing. The index was calculated as 
	the sum of the five items (range 5–18). It was computed only for respondents with 
	complete data on all five items (n = 34,533; 93.6% of the sample).

	An initial index included belief in God, change in spirituality over the life course, 
	and the view that science and religion are compatible. That version had lower 
	reliability (α = .73). The life-course item also could not be scored for the 24% of 
	respondents who reported becoming both more and less spiritual at different times. 
	Belief in God was excluded because it forms part of the religious belief index. Two 
	experience items were added because their content refers directly to spiritual peace 
	and to a presence beyond this world.

	The final index demonstrated acceptable internal consistency (Cronbach's α = .75; 
	Nunnally, 1978). Corrected item–total correlations ranged from .48 (belief in a soul) 
	to .67 (feeling the presence of something beyond this world). Removing belief in a soul 
	would have raised alpha only marginally (.76), so the item was retained for its 
	theoretical relevance. Scores were negatively skewed (M = 13.62, SD = 3.59, Mdn = 14, 
	skewness = −0.55). At least 10% of respondents obtained the maximum score, indicating a 
	ceiling. Reliability and descriptive statistics are unweighted. */
	
		

* Items in the current indices that are not belief items - can serve as validation 	variables after the LPA. For example, if the Spiritual but not religious profile is valid, its members should be more likely than other profiles to report becoming less religious and more spiritual. For that use, keep them as categorical variables with all four responses, including "some of both". 

	codebook lifespir lifedir
	
	

**# 3. Religious tradition - current and raised

	tab1 reltrad frmrel
	
	recode reltrad (1100 1200 1300=1 "Protestant") (10000=2 "Catholic") ///
		(20000 30000 40001 40002=3 "Other Christian") (50000=4 "Jewish") ///
		(60000=5 "Muslim") (70000=6 "Buddhist") (80000=7 "Hindu") ///
		(90001 90002=8 "Other religions") ///
		(100000=9 "Religiously unaffiliated"), gen(rel_curr9)
	la var rel_curr9 "Current religious tradition (9 categories)"

	recode frmrel (1000=1 "Protestant") (10000=2 "Catholic") ///
		(20000 30000 40001 40002=3 "Other Christian") (50000=4 "Jewish") ///
		(60000=5 "Muslim") (70000=6 "Buddhist") (80000=7 "Hindu") ///
		(90001 90002=8 "Other religions") ///
		(100000=9 "Religiously unaffiliated"), gen(rel_raised9)
	la var rel_raised9 "Religious tradition raised in (9 categories)"

	tab reltrad rel_curr9, mi
	tab frmrel rel_raised9, mi



**# 4. Demographics

	codebook birthdecade gender educrec inc_sdt1 ideo party partyln usgen marital

	// Age group (reversed birth decade: 1 = youngest, 7 = oldest)

	recode birthdecade (7=1 "Born 2000s (about 18-24)") ///
	(6=2 "Born 1990s (about 24-34)") ///
    (5=3 "Born 1980s (about 34-44)") (4=4 "Born 1970s (about 44-54)") ///
    (3=5 "Born 1960s (about 54-64)") (2=6 "Born 1950s (about 64-74)") ///
    (1=7 "Born 1940s or earlier (about 74+)"), gen(agegrp)
	la var agegrp "Age group (from birth decade; higher = older)"

	* Gender
	recode gender (1=0 "Man") (2=1 "Woman") (3=.), gen(woman)
	la var woman "Gender (0 = man, 1 = woman)"

	* Education: no recode needed
	la var educrec "Education (4 categories)"

	* Family income, 4 categories
	recode inc_sdt1 (1/2=1 "Less than $40,000") (3/4=2 "$40,000 to less than $75,000") ///
		(5/6=3 "$75,000 to less than $125,000") (7/8=4 "$125,000 or more"), gen(income4)
	la var income4 "Family income (4 categories)"

	* Ideology: no recode needed (1 = very conservative, 5 = very liberal)
	la var ideo "Political ideology (1 = very conservative, 5 = very liberal)"

	* Party identification including leaners, 3 categories
	gen partyid3 = .
	replace partyid3 = 1 if party == 1 | partyln == 1
	replace partyid3 = 2 if party == 2 | partyln == 2
	replace partyid3 = 3 if inlist(party, 3, 4) & missing(partyln)
	la def partyid3 1 "Republican or lean Republican" 2 "Democrat or lean Democrat" ///
		3 "Independent or other, no lean"
	la val partyid3 partyid3
	la var partyid3 "Party identification including leaners"

	* Immigrant generation: no recode needed
	la var usgen "Immigrant generation (3 categories)"

	* Marital status, 4 categories
	recode marital (1 2=1 "Married or living with a partner") ///
		(3 4=2 "Divorced or separated") ///
		(5=3 "Widowed") (6=4 "Never married"), gen(marital4)
	la var marital4 "Marital status (4 categories)"

	* Checks
	tab birthdecade agegrp, mi
	tab gender woman, mi
	tab educrec, mi
	tab inc_sdt1 income4, mi
	tab ideo, mi
	tab party partyln, mi
	tab partyid3, mi
	tab usgen, mi
	tab marital marital4, mi

	
	* Race and ethnicity (racecmb, hisp) are not used in this analysis.


**# 5. Childhood socialization

	codebook chimprel_a chimprel_b chattend chedu1 chedu2


	// Importance of religion when growing up (reverse: 0 = not at all, 3 = very important)
	recode chimprel_a (4=0 "Not at all important") (3=1 "Not too important") ///
    (2=2 "Somewhat important") (1=3 "Very important"), gen(ch_imp_self)
	la var ch_imp_self "Importance of religion to self when growing up (0-3)"

	recode chimprel_b (4=0 "Not at all important") (3=1 "Not too important") ///
    (2=2 "Somewhat important") (1=3 "Very important"), gen(ch_imp_family)
	la var ch_imp_family "Importance of religion to family when growing up (0-3)"

	* Childhood religiosity scale (mean of the two items, complete cases only)
	egen ch_nmiss = rowmiss(ch_imp_self ch_imp_family)
	egen child_relig = rowmean(ch_imp_self ch_imp_family) if ch_nmiss == 0
	la var child_relig "Childhood religiosity (0-3, mean of 2 items)"

	* Childhood attendance (reverse: 0 = never, 5 = more than once a week)
	recode chattend (6=0 "Never") (5=1 "Seldom") (4=2 "A few times a year") ///
		(3=3 "Once or twice a month") (2=4 "Once a week") ///
		(1=5 "More than once a week"), gen(ch_attend)
	la var ch_attend "Attended religious services when growing up (0-5)"

	* Religious schooling (0 = none)
	recode chedu1 (1=0 "None") (2=1 "1-3 years") (3=2 "4-6 years") ///
	    (4=3 "7 years or more"), gen(ch_relschool)
	la var ch_relschool "Years at a private religious school (K-12)"

	recode chedu2 (1=0 "None") (2=1 "1-3 years") (3=2 "4-6 years") ///
		(4=3 "7 years or more"), gen(ch_religed)
	la var ch_religed "Years of Sunday school or other religious education"

	* Checks
	tab chimprel_a ch_imp_self, mi
	tab chimprel_b ch_imp_family, mi
	tab chattend ch_attend, mi
	tab chedu1 ch_relschool, mi
	tab chedu2 ch_religed, mi

	* Diagnostics
	tab ch_nmiss
	alpha ch_imp_self ch_imp_family
	summarize child_relig, detail
	pwcorr child_relig ch_attend ch_relschool ch_religed, obs
		
	
	
**# 6. Adult relational context

	codebook relmatch sprelimp 
	
	
	// Religious match with spouse or partner
	gen partner_rel = .
	replace partner_rel = 0 if inlist(marital, 3, 4, 5, 6)
	replace partner_rel = 1 if inlist(marital, 1, 2) & relmatch == 0
	replace partner_rel = 2 if inlist(marital, 1, 2) & relmatch == 1
	la def partner_rel 0 "No spouse or partner" 1 "Spouse/partner same religion" ///
		2 "Spouse/partner different religion"
	la val partner_rel partner_rel
	la var partner_rel "Religious match with spouse or partner"

	// Importance of religion to spouse or partner
	gen sp_relimp = .
	replace sp_relimp = 0 if inlist(marital, 3, 4, 5, 6)
	replace sp_relimp = 1 if inlist(marital, 1, 2) & sprelimp == 4
	replace sp_relimp = 2 if inlist(marital, 1, 2) & sprelimp == 3
	replace sp_relimp = 3 if inlist(marital, 1, 2) & sprelimp == 2
	replace sp_relimp = 4 if inlist(marital, 1, 2) & sprelimp == 1
	la def sp_relimp 0 "No spouse or partner" 1 "Not at all important" ///
		2 "Not too important" 3 "Somewhat important" 4 "Very important"
	la val sp_relimp sp_relimp
	la var sp_relimp "Importance of religion to spouse or partner"

	// Checks
	tab marital partner_rel, mi
	tab relmatch partner_rel, mi
	tab marital sp_relimp, mi
	tab sprelimp sp_relimp, mi

	
	
	
**# 7. Wellbeing 

	codebook happy satis_a satis_b


	* i. Reverse-code so 0 = lowest, higher = better
	recode happy (3=0 "Not too happy") (2=1 "Pretty happy") ///
		(1=2 "Very happy"), gen(wb_happy)
	la var wb_happy "General happiness (0-2)"

	recode satis_a (5=0 "Poor") (4=1 "Fair") (3=2 "Good") (2=3 "Very good") ///
		(1=4 "Excellent"), gen(wb_health)
	la var wb_health "Self-rated health (0-4)"

	recode satis_b (5=0 "Poor") (4=1 "Fair") (3=2 "Good") (2=3 "Very good") ///
		(1=4 "Excellent"), gen(wb_family)
	la var wb_family "Quality of family life (0-4)"

	* Checks
	tab happy wb_happy, mi
	tab satis_a wb_health, mi
	tab satis_b wb_family, mi

	* Index (raw sum, complete cases only)
	local wbitems wb_happy wb_health wb_family
	egen wb_nmiss = rowmiss(`wbitems')
	egen wellbeing = rowtotal(`wbitems') if wb_nmiss == 0
	la var wellbeing "Wellbeing index (sum of 3 items, 0-10)"

	* Diagnostics
	tab wb_nmiss
	alpha `wbitems', item
	summarize wellbeing, detail
	
	


	**# 8. Analytic sample flag and random halves (factor analysis is done in R)

	egen lpa_nmiss = rowmiss(religbehav religbelif spirbehav spirbelif)
	gen lpa_sample = (lpa_nmiss == 0)
	tab lpa_sample

	sort respid
	set seed 20261016
	gen half = cond(runiform() < .5, 1, 2) if lpa_sample == 1
	tab half
	
	
	
	log close
	

		
save "path/to/data/RLS Public Use File_2_Clean.dta", replace


	
	