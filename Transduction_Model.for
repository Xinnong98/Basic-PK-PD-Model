**********************************************************************
C                           ADAPT                                     *
C                         Version 5                                   *
C**********************************************************************
C                                                                     *
C                           MODEL                                     *
C                                                                     *
C    This file contains Fortran subroutines into which the user       *
C    must enter the relevant model equations and constants.           *
C    Consult the User's Guide for details concerning the format for   *
C    entered equations and definition of symbols.                     *
C                                                                     *
C       1. Symbol-  Parameter symbols and model constants             *
C       2. DiffEq-  System differential equations                     *
C       3. Output-  System output equations                           *
C       4. Varmod-  Error variance model equations                    *
C       5. Covmod-  Covariate model equations (ITS,MLEM)              *
C       6. Popinit- Population parameter initial values (ITS,MLEM)    *
C       7. Prior -  Parameter mean and covariance values (ID,NPD,STS) *
C       8. Sparam-  Secondary parameters                              *
C       9. Amat  -  System state matrix                               *
C                                                                     *
C**********************************************************************

C######################################################################C

        Subroutine SYMBOL
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

CC
C----------------------------------------------------------------------C
C   Enter as Indicated                                                 C
C----c-----------------------------------------------------------------C

      NDEqs   =  18   ! Enter # of Diff. Eqs.
      NSParam =  10   ! Enter # of System Parameters.
      NVparam =  2   ! Enter # of Variance Parameters.
      NSecPar =  0   ! Enter # of Secondary Parameters.
      NSecOut =  0  ! Enter # of Secondary Outputs (not used).
      Ieqsol  =  1  ! Model type: 1 - DIFFEQ, 2 - AMAT, 3 - OUTPUT only.
      Descr   = ' Indirect Binding  '

CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each System Parameter (eg. Psym(1)='Kel')         C
C----c-----------------------------------------------------------------C

		Psym(1)='cl'
		Psym(2)='cld'
		Psym(3)='v2'
		Psym(4)='V1'
		Psym(5)='ti'
		Psym(6)='tau'
		Psym(7)='Em'
		Psym(8)='EC50'
          Psym(9)='ga'
		  Psym(10)='bl'

		  

		


CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each Variance Parameter {eg: PVsym(1)='Sigma'}    C
C----c-----------------------------------------------------------------C
		PVsym(1)='intercept'
		PVsym(2)='Sigma'
		
CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each Secondary Parameter {eg: PSsym(1)='CLt'}     C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine DIFFEQ(T,X,XP)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 T,X(MaxNDE),XP(MaxNDE)
		Real*8 cl,cld,v2,V1,ti,Em,EC50,ga,tau
C----------------------------------------------------------------------C
C   Enter Differential Equations Below  {e.g.  XP(1) = -P(1)*X(1) }    C
C----c-----------------------------------------------------------------C
		cl=P(1)
		cld=P(2)
		v2=P(3)
		V1=P(4)
		ti=P(5)
		tau=P(6)
		Em=P(7)
		EC50=P(8)
          ga=P(9)

		if (T.le.ti)then
		XP(1) = -10000000/ti
		XP(2) = 10000000/ti-(cl/v1+cld/v1)*X(2)+(cld/v2*X(3))
		XP(3) = cld/v1*X(2)-(cld/v2*X(3))
		
		XP(4) = -20000000/ti
		XP(5) = 20000000/ti-(cl/v1+cld/v1)*X(5)+(cld/v2*X(6))
		XP(6) = cld/v1*X(5)-(cld/v2*X(6))
		
		XP(7) = -30000000/ti
		XP(8) = 30000000/ti-(cl/v1+cld/v1)*X(8)+(cld/v2*(9))
		XP(9) = cld/v1*X(8)-(cld/v2*X(9))
		else
		XP(1) = 0
		XP(2) = -(cl/v1+cld/v1)*X(2)+(cld/v2*X(3))
		XP(3) = cld/v1*X(2)-(cld/v2*X(3))
		
		XP(4) = 0
		XP(5) = -(cl/v1+cld/v1)*X(5)+(cld/v2*X(6))
		XP(6) = cld/v1*X(5)-(cld/v2*X(6))
		
		XP(7) = 0
		XP(8) = -(cl/v1+cld/v1)*X(8)+(cld/v2*X(9))
		XP(9) = cld/v1*X(8)-(cld/v2*X(9))
		end if
		XP(10)=1/tau*((Em*(X(2)/V1)**ga/(EC50**ga+(X(2)/V1)**ga))-X(10))
		XP(11)=1/tau*(X(10)-X(11))
		XP(12)=1/tau*(X(11)-X(12))
		XP(13)=1/tau*((Em*(X(5)/V1)**ga/(EC50**ga+(X(5)/V1)**ga))-X(13))
		XP(14)=1/tau*(X(13)-X(14))
		XP(15)=1/tau*(X(14)-X(15))
		XP(16)=1/tau*((Em*(X(8)/V1)**ga/(EC50**ga+(X(8)/V1)**ga))-X(16))
		XP(17)=1/tau*(X(16)-X(17))
		XP(18)=1/tau*(X(17)-X(18))
		
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine OUTPUT(Y,T,X)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 Y(MaxNOE),T,X(MaxNDE),V1,bl

CC
C----------------------------------------------------------------------C
C   Enter Output Equations Below   {e.g.  Y(1) = X(1)/P(2) }           C
C----c-----------------------------------------------------------------C
		V1=P(4)
		bl=P(10)
        Y(1) = X(2)/V1
		Y(2) = X(5)/V1
		Y(3) = X(8)/V1
		Y(4) = bl-X(12)
		Y(5) = bl-X(15)
		Y(6) = bl-X(18)
		Y(7) = bl


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine VARMOD(V,T,X,Y)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 V(MaxNOE),T,X(MaxNDE),Y(MaxNOE)

CC
C----------------------------------------------------------------------C
C   Enter Variance Model Equations Below                               C
C         {e.g. V(1) = (PV(1) + PV(2)*Y(1))**2 }                       C
C----c-----------------------------------------------------------------C
		 V(1) = (PV(1) + PV(2)*Y(1))**2 
		 V(2) = (PV(1) + PV(2)*Y(2))**2
		 V(3) = (PV(1) + PV(2)*Y(3))**2 
		 V(4) = (PV(1) + PV(2)*Y(4))**2 
		 V(5) = (PV(1) + PV(2)*Y(5))**2
		 V(6) = (PV(1) + PV(2)*Y(6))**2
		 V(7) = (PV(1) + PV(2)*Y(7))**2


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine COVMOD(Pmean, ICmean, PC)
C  Defines any covariate model equations (MLEM, ITS)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 PC(MaxNCP)
        Real*8 Pmean(MaxNSP+MaxNDE), ICmean(MaxNDE)

CC
C----------------------------------------------------------------------C
C     Enter # of Covariate Parameters                                  C
C----c-----------------------------------------------------------------C

        NCparam = 0    ! Enter # of Covariate Parameters.

CC
C----------------------------------------------------------------------C
C   Enter Symbol for Covariate Params {eg: PCsym(1)='CLRenal'}         C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   For the Model Params. that Depend on Covariates Enter the Equation C
C         {e.g. Pmean(1) =  PC(1)*R(2) }                               C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine POPINIT(PmeanI,ICmeanI,PcovI,ICcovI, PCI)
C  Initial parameter values for population program parameters (ITS, MLEM)

        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 PmeanI(MaxNSP+MaxNDE), ICmeanI(MaxNDE)
        Real*8 PcovI(MaxNSP+MaxNDE,MaxNSP+MaxNDE), ICcovI(MaxNDE,MaxNDE)
        Real*8 PCI(MaxNCP)

CC
C----------------------------------------------------------------------C
C  Enter Initial Values for Population Means                           C
C          {  e.g. PmeanI(1) = 10.0    }                               C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Initial Values for Pop. Covariance Matrix (Lower Triang.)    C
C         {  e.g. PcovI(2,1) = 0.25    }                               C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Values for Covariate Model Parameters                        C
C         {  e.g. PCI(1) = 2.0    }                                    C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine PRIOR(Pmean,Pcov,ICmean,ICcov)
C  Parameter mean and covariance values for MAP estimation (ID,NPD,STS)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 Pmean(MaxNSP+MaxNDE), ICmean(MaxNDE)
        Real*8 Pcov(MaxNSP+MaxNDE,MaxNSP+MaxNDE), ICcov(MaxNDE,MaxNDE)

CC
C----------------------------------------------------------------------C
C  Enter Nonzero Elements of Prior Mean Vector                         C
C          {  e.g. Pmean(1) = 10.0    }                                C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Nonzero Elements of Covariance Matrix (Lower Triang.)       C
C         {  e.g. Pcov(2,1) = 0.25    }                                C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine SPARAM(PS,P,IC)
        Implicit None

        Include 'globals.inc'

        Real*8 PS(MaxNSECP), P(MaxNSP+MaxNDE), IC(MaxNDE) 

CC
C----------------------------------------------------------------------C
C   Enter Equations Defining Secondary Paramters                       C
C           {  e.g.  PS(1) = P(1)*P(2)   }                             C
C----c-----------------------------------------------------------------C

        
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End
	  	  
C######################################################################C

        Subroutine AMAT(A)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 A(MaxNDE,MaxNDE)

        DO I=1,Ndeqs
           Do J=1,Ndeqs
              A(I,J)=0.0D0
           End Do
        End Do

CC
C----------------------------------------------------------------------C
C   Enter non zero elements of state matrix  {e.g.  A(1,1) = -P(1) }   C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C