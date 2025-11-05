================================================================
=== readme.txt =================================================
================================================================

Provided are as part of the code release the configuration files necessary to run the key model experiments presented in the paper.
The intention is to provide an oppertunity to question the paper assumptions and interpretation through re-analysis,
as well as the creation of new and different experiments. (Plus, to provide a means to replicate published results.)
This readme file details how the experiments can be run.
Refer to the cookie manual:
https://github.com/genie-model/cookiedoc
for details on model code installation and configuration, locating and visualizing model results, etc.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PUBLICATION DETAILS
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Inclusion of MyAMI-derived Mg/Ca corrections to the marine carbonate system in the cGENIE.cookie Earth system model (v.0.9) 

Markus Adloff, Terra Ganey, Mathis P. Hain, Michael J. Henehan, Sarah E. Greene, Andy Ridgwell


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
EDITING LOG
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

2025/11/04 -- README.txt file created by [AR]

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUMMARY OF EXPERIMENTS
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The names of all experiments start:
cookie.CBSR.p_worbe2.BASES
and then are followed by 3 (or 4) underscore-seperated codes -- 'ID' in the table below.
Note that not all experiments are used in the paper. The relevant Figure(s) is indicated.
TS == truncated T,S limits

========================================================================================================
Ca/Mg	CORRECTION	SYSTEM CONFIG	SEDIMENTS	WEATHERING	ID				Table 1 		FIGURE(S)
========================================================================================================
modern	TZ			ocean-only		no			none		m-TZ-ocean		PI-A-ocean		3,SI.1,SI.2
		MY			ocean-only		no			none		m-MY-ocean		PI-B-ocean		3,SI.1
		NONE		ocean-only		no			none		m-NO-ocean		(not used)
--------------------------------------------------------------------------------------------------------
Eocene	TZ			ocean-only		no			none		e-TZ-ocean		(not used)
		MY			ocean-only		no			none		e-MY-ocean		EE-B-ocean		2
		NONE		ocean-only		no			none		e-NO-ocean		(not used)
========================================================================================================
modern	TZ			closed			yes			tracking	m-TZ-closed		PI-A-closed		5,7
		MY			closed			yes			tracking	m-MY-closed		PI-B-closed		5,7
		NONE		closed			yes			tracking	m-NO-closed		PI-C-closed		5,7
--------------------------------------------------------------------------------------------------------
Eocene	TZ			closed			yes			tracking	e-TZ-closed		EE-A-closed		5,7
		MY			closed			yes			tracking	e-MY-closed		EE-B-closed		5,6,7
		NONE		closed			yes			tracking	e-NO-closed		EE-C-closed		5,6,7
========================================================================================================
modern	TZ			open			yes			fixed		m-TZ-open		PI-A-open		4,7
		MY			open			yes			fixed		m-MY-open		PI-B-open		4,7
		NONE		open			yes			fixed		m-NO-open		PI-C-open		7
--------------------------------------------------------------------------------------------------------
Eocene	TZ			open			yes			fixed		e-TZ-open		EE-A-open		7,8,SI.3,SI.4,SI.5
		MY			open			yes			fixed		e-MY-open		EE-B-open		7,8,SI.3,SI.4,SI.5,SI.6
		NONE		open			yes			fixed		e-NO-open		EE-C-open		7
========================================================================================================
modern	MY-TS		ocean-only		no			none		m-MY-TS-ocean	PI-B2-ocean		SI.2
Eocene	MY-TS		ocean-only		no			none		e-MY-TS-ocean	(not used)
========================================================================================================

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RUNNING THE EXPERIMENTS
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In the following example commands for running the model:
x 	== 	e, m
YY 	== 	NO, TZ, MY

(1)	For all the 'ocean-only' marine carbon cycle experiments -- to run a 10,000 year spin-up:

./runcookie.sh cookie.CB.p_worbe2.BASES PUBS\submitted\Adloff_etal.2025 cookie.CB.p_worbe2.BASES.x_YY_ocean.SPIN 10000

(2)	For the 'closed system' ocean+sediment carbon cycle configuration -- to run a 20,000 year spin-up:

./runcookie.sh cookie.CBSR.p_worbe2.BASES PUBS\submitted\Adloff_etal.2025 cookie.CBSR.p_worbe2.BASES.x_YY_closed.SPIN 20000

(2)	For the 'open system' ocean+sediment carbon cycle configuration -- to run a 100,000 year spin-up:

./runcookie.sh cookie.CBSR.p_worbe2.BASES PUBS\submitted\Adloff_etal.2025 cookie.CBSR.p_worbe2.BASES.x_YY_open.SPIN 100000

(4)	Finally, for the MyAmi Mg/Ca correction (ocean-only configuration) but with default T and S limits imposed:

./runcookie.sh cookie.CB.p_worbe2.BASES PUBS\submitted\Adloff_etal.2025 cookie.CB.p_worbe2.BASES.x_MY_TS_ocean.SPIN 10000

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Note that all experiments are run from:
$HOME/cgenie.cookie/genie-main
(unless a different installation directory has been used)

================================================================
================================================================
================================================================
