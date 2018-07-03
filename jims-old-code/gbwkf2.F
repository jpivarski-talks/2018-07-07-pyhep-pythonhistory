	DOUBLE PRECISION FUNCTION GBWKF(RM,GAM,WSPREAD,HC)
        IMPLICIT NONE
c 	implicit real (A-H,O-Z), integer (I-N)

C   Observed shape of resonance peak, starting with Breit-Wigner
C   Normalized per unit W, MeV^-1
C   F_KF(x,s) convoluted with unit-area Breit-Wigner, and with beam resolution
C    unit-area Gaussian
C   If GAM<1 keV, the Breit-Wigner is replaced with a delta function
C   BW(w) = (GAM/2*pi)/[(w-M)^2+GAM^2/4] (unit integral over w)
C   Gaussian(w) = exp{-.5*[(w-M)^2/WSPREAD^2]}/(2*pi*WSPREAD) (unit int...)
C   Calls FUNCTION FKF to get Kuraev-Fadin radiative tail function, F_KF(x,s)
C   BW*KF -> BWKF; G*BWKF -> GBWKF
C   All energies are in MeV
C   RM = central resonance mass
C   GAM= resonance Breit-Wigner total width
C   WSPREAD = r.m.s. c.m. beam energy spread, in MeV
C   WT = true c.m. e+e- collision energy for particular event
C   WR = rest energy of produced resonance for particular event
C   WA = WT*RM/WR
C   WC = centroid of CESR energy distribution
C   HC = WC - RM  (and so on for other H...)
C
        INTEGER   MKBIN, MTBIN, JTBIN
        REAL*8    ROOT2, PI, RTWOPI, RTWOBPI, PMAX, TINY

	PARAMETER (MKBIN=1500)                ! # bins in WA-RM
	PARAMETER (MTBIN=2000)                ! # bins of HT=WT-RM, each sign
	PARAMETER (ROOT2=1.414213d0)
	PARAMETER (PI=3.1415926d0)
	PARAMETER (RTWOPI=2.506628d0)           ! sqrt(2pi)
	PARAMETER (RTWOBPI=0.7978845d0)         ! sqrt(2/pi)
	PARAMETER (JTBIN=20)                  ! # HT bins within GAM
	PARAMETER (PMAX=8.d0)                   ! max power in Gaussian
	PARAMETER (TINY=1.d-3)                ! minimum GAM used
	REAL*8    FKFT(MKBIN)	              ! F_KF(x,s) (Kuraev-Fadin eq.28)
	REAL*8    HA(MKBIN)                   ! WA-RM of centroid WA bin
	REAL*8    DWA(MKBIN)                  ! bin width in WA-RM
	REAL*8    HT(0:MTBIN)
	REAL*8    DWT(0:MTBIN)                ! WT-RM bin center
	REAL*8    BWKF(-MTBIN:MTBIN)          ! result of BW-F_KF convolution
        REAL*8    RM, GAM, WSPREAD, HC
        REAL*8    RMSAV, GAMSAV, DWA1, HA1
        INTEGER   KBIN, J, ITBIN
        REAL*8    X, S
        REAL*8    AJTB, DWT1, HT1, D
        REAL*8    GAMBPIM, GAMSB4, H, WT, SUM, SIGS, P, WA, HR
        REAL*8    FKF
        EXTERNAL  FKF
	DATA RMSAV,GAMSAV/0.d0,0.d0/
cww     DATA WSSAV/0./
C
C---Is this a new set of parameters (RM,GAM,WSPREAD)?
cww	IF(WSPREAD.NE.WSSAV)GO TO 18

	IF(RM.NE.RMSAV)GO TO 18
	IF(GAM.NE.GAMSAV)GO TO 28  !JRP!!! Used to be GOTO 18
	GO TO 19
C
C---Save new parameters
 18	GAMSAV=GAM
	RMSAV=RM
cww	WSSAV=WSPREAD
C	
C---Initialize table of F_KF values
C     F_KF is calculated once at an array of energies WA, and stored as FKFT.
C     F_KF is a function of x, and only weakly (in logs) a function of s.
C     Nominally, s = WT^2 and x = (s-WR^2)/s.
C     We approximate F_KF as a function of a single variable WA = sqrt(s').
C     We take s' = (WT*RM/WR)^2 and x' = (s'-RM^2)/s', which differs only
C      in the explicit dependence of F_KF on s in log terms.
C     
	DWA1=.5/MKBIN                       ! width of bin1=(1 MeV)/MKBIN
	DWA(1)=DWA1
	HA(1)=.5d0*DWA1                       ! middle of bin1
	HA1=0.                              ! bottom of bin1
	DO 10 KBIN=1,MKBIN                  ! loop over WA's
	  IF(KBIN.GT.1)THEN
	    J=KBIN-1
	    HA1=HA1+DWA(J)                  ! bottom of bin
	    DWA(KBIN)=KBIN*DWA1             ! bin width (prop to WA-RM)
cw	    HA(KBIN)=DWA(KBIN)/DLOG(1.+DWA(KBIN)/HA1)   ! 'centroid' of bin
	    HA(KBIN)=HA1+.5*DWA(KBIN)       ! middle of bin
	  ENDIF
 	  S=(HA(KBIN)+RM)**2
	  X=(S-RM**2)/S
	  FKFT(KBIN)=FKF(X,S)               ! store F_KF value
 10	CONTINUE
ct	type 1000,MKBIN,HA(MKBIN),FKFT(1)   ! diagnostic
ct 1000	format(' MKBIN,HA(MKBIN),FKFT(1)=',I6,F12.1,E12.4)
C
C---Define bins in HT
 28	CONTINUE
	IF(GAM.LT.TINY)GO TO 20             ! skip BW convolution?
	AJTB=FLOAT(JTBIN)
	DWT1=GAM/AJTB                       ! inside binwid in MeV
	DWT(0)=DWT1                         ! width of bin0
	HT1=-.5d0*DWT1                        ! bottom of bin0
	HT(0)=0.d0                            ! middle of bin0
	DO 17 ITBIN=1,MTBIN                 ! loop over WT's
	  IF(ITBIN.LE.JTBIN)THEN            ! inside GAM?
	    D=DWT1                          ! width of bin (fixed)
	    HT1=HT1+DWT1                    ! bottom of bin
	    HT(ITBIN)=HT1+.5d0*DWT1           ! middle of bin
	  ELSE                              ! outside GAM
	    D=DWT1*FLOAT(ITBIN)/AJTB        ! width of bin (prop to WT-RM)
	    HT1=HT1+DWT(ITBIN-1)            ! bottom of bin
cw	    HT(ITBIN)=SQRT(D/(1./HT1-1./(HT1+D)))  ! 'centroid' of bin
	    HT(ITBIN)=HT1+.5d0*DWT(ITBIN)     ! middle of bin
	  ENDIF
	  DWT(ITBIN)=D
 17	CONTINUE
C
C---Load table of BW convoluted with F_KF, for an array of WT values
	GAMBPIM=GAM/(PI*RM)                ! GAM/(pi*RM)
	GAMSB4=GAM**2/4.d0
	DO 14 ITBIN=-MTBIN,MTBIN           ! loop over WT's
	  H=HT(IABS(ITBIN))                ! middle of bin
	  IF(ITBIN.LT.0)H=-H
	  WT=H+RM
	  SUM=0.d0
	  DO 15 KBIN=1,MKBIN               ! sum over WA's
	    WA=HA(KBIN)+RM
	    HR=(WT/WA-1.)*RM
cg	    SUM=SUM+FKFT(KBIN)*DWA(KBIN)/(HR**2+GAMSB4)  ! BW factor
	    SUM=SUM+FKFT(KBIN)*DWA(KBIN)/(HR**2+GAMSB4*(1.+HR/RM)**2)
 15	  CONTINUE
	  IF(DABS(HR)/GAM.GT.8.d0)THEN        ! have we got all of BW?
	    BWKF(ITBIN)=GAMBPIM*SUM        ! BW*KF convolution
C            2/RM converts variable of integration from x to W.
C            GAM/(2pi) normalizes Breit-Wigner
	  ELSE
	    S=WT**2
	    BWKF(ITBIN)=2.d0*FKF((S-RM**2)/S,S)/RM  ! approx BW as delta fn
	  ENDIF
 14	CONTINUE
ct	type 1003,mtbin,ht(mtbin),bwkf(0)  ! diagnostic
ct 1003	format(' MTBIN,HT(MTBIN),BWKF(0)=',I6,F12.1,E12.4)
C
C---Calculate [BW convoluted with F_KF] convoluted with Gaussian, at given HC
 19	IF(GAM.LT.TINY)GO TO 20            ! skip BW convolution?
	SUM=0.d0
	DO 11 ITBIN=-MTBIN,MTBIN           ! sum over WT's
	  H=HT(IABS(ITBIN))                ! WT-RM
	  IF(ITBIN.LT.0)H=-H
	  SIGS=(H-HC)/(ROOT2*WSPREAD)      ! sqrt of exponent
	  IF((ITBIN.EQ.MTBIN).AND.(SIGS.LT.2.))GO TO 12 ! running out of bins?
	  IF(SIGS.GT.3.)GO TO 13           ! Gaussian getting negligible?
	  P=SIGS**2                        ! exponent
	  IF(P.LT.PMAX)THEN
	    IF(DWT(IABS(ITBIN)).GT..2*WSPREAD)GO TO 12  ! binning too coarse?
            SUM=SUM+EXP(-P)*BWKF(ITBIN)*DWT(IABS(ITBIN))
	  ENDIF
 11	CONTINUE
 13	GBWKF=SUM/(RTWOPI*WSPREAD)         ! convolution
C         1/(sqrt(2*pi)*WSPREAD) normalizes Gaussian
	RETURN
C
C---Don't smear if overflow MKBINs
 12	S=(HC+RM)**2
	GBWKF=2.d0*FKF((S-RM**2)/S,S)/RM     ! ignore GAM, WSPREAD
	RETURN
C
C---Gaussian, KF convolution in case of negligible GAM
 20	SUM=0.
	DO 21 KBIN=1,MKBIN                     ! loop over array points
	  SIGS=(HA(KBIN)-HC)/(ROOT2*WSPREAD)
	  IF((KBIN.EQ.MKBIN).AND.(SIGS.LT.2.))GO TO 12 ! running out of bins?
	  IF(SIGS.GT.3.)GO TO 23              ! Gaussian getting negligible?
	  P=SIGS**2
	  IF(P.LT.PMAX)THEN
	    IF(DWT(IABS(ITBIN)).GT..2*WSPREAD)GO TO 12  ! binning too coarse?
	    SUM=SUM+EXP(-P)*FKFT(KBIN)*DWA(KBIN)
	  ENDIF
 21	CONTINUE
 23	GBWKF=RTWOBPI*SUM/(RM*WSPREAD)
C         2/RM converts variable of integration from x to W.
C         1/(sqrt(2*pi)*WSPREAD) normalizes Gaussian
	RETURN
	END
C
C--------------------------------------------------------------------------
	DOUBLE PRECISION FUNCTION FKF(X,S)
        IMPLICIT NONE
C   Radiative tail spectrum from Kuraev and Fadin (equation 28)
C   Does not include contribution of radiated e+e- pairs (theta term)
C    (Use FUNCTION FKFP to get same F with pairs included)
C   This routine assumes that resonance Gamma_ee ~ 0
C   All energies are in MeV
C   RM = central resonance mass
C   WT = true c.m. e+e- energy for particular event
C   WR = rest energy of produced resonance for particular event
C    assumed = RM for this routine
C   WC = centroid of CESR energy spectrum
C   S = WT^2
C   X = (WT^2 - WR^2)/WT^2   (~DeltaE/E)
C   FKF diverges at x->0 like 1/(x^.9), so we average over x for x<XAV
C   DO NOT call this function for x<=0; it's not protected!
        REAL*8 XAV, XFAST, ALPBPI, PI, AM2, FAC
	PARAMETER (XAV=1.d-6)                   ! average over x<XAV
	PARAMETER (XFAST=.01d0)                 ! approx logs for x<XFAST
	PARAMETER (ALPBPI=2.322819d-3)		! alpha/pi
	PARAMETER (PI=3.1415926d0)		! pi, pi^2
	PARAMETER (AM2=.2611199d0)		! m_e^2 in MeV^2
	PARAMETER (FAC=2.78987d0)               ! pi^2/3 - .5
        REAL*8 X, S

        REAL*8 ALOGRF, ALOGR, BETAF, BETA, BETA2, BETA2F
        REAL*8 F1, FKFAV, T1, T2, T3
	LOGICAL LFIRST, LFAST
	DATA LFIRST/.TRUE./                     ! T: haven't hit singularity
	DATA LFAST/.TRUE./                      ! T: haven't done low-x logs
C
C---Initialize logs for fast calculation at low x
	IF(X.LE.XFAST)THEN                      ! x low enough to approx logs?
	  IF(LFAST)THEN                         ! no low-x logs yet?
	    ALOGRF=DLOG(S/AM2)			! "L"  (K&F's notation)
	    BETAF=2.d0*ALPBPI*(ALOGRF-1.D0)	! "beta" (same as J&S's "t")
	    BETA2F=BETAF**2
	    LFAST=.FALSE.
	  ENDIF
	  ALOGR=ALOGRF                          ! use approximate low-x logs
	  BETA=BETAF
	  BETA2=BETA2F
	ELSE
	  ALOGR=DLOG(S/AM2)			! slow calculation for logs
	  BETA=2.d0*ALPBPI*(ALOGR-1.)
	  BETA2=BETA**2
	ENDIF
C
C---Average over small-x singularity
	IF(X.LE.XAV)THEN                        ! is x at the singularity?
	  IF(LFIRST)THEN                        ! first time at singularity?
	    F1=1.d0+.75d0*BETA+ALPBPI*FAC       ! integrate over singularity
     2       +BETA**2*(-.4370503-ALOGR/72.)
	    FKFAV=F1*XAV**(BETA-1.)
	    LFIRST=.FALSE.
	  ENDIF
	  FKF=FKFAV                             ! FKF = saved singularity ave.
	  RETURN
	ENDIF
C
C---For non-infinitesimal x
	T1=BETA*X**(BETA-1.d0)*(1.d0+ALPBPI*FAC+.75d0*BETA ! 1st term in eq. 28
     2   -(BETA2/24.d0)*(ALOGR/3.d0+10.4892d0))
	T2=BETA*(1.d0-.5d0*X)                                ! 2nd term ..
	T3=.125*BETA2*
     2    (-4.d0*(2.d0-X)*DLOG(X)-((1.d0+3.d0*(1.d0-X**2))/X)  ! 3rd term
     2   *DLOG(1.d0-X)-6.d0+X)
	FKF=T1-T2+T3
	RETURN
	END
