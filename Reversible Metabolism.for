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

      NDEqs   =  8   ! Enter # of Diff. Eqs.
      NSParam =  10   ! Enter # of System Parameters.
      NVparam =  2   ! Enter # of Variance Parameters.
      NSecPar =  0   ! Enter # of Secondary Parameters.
      NSecOut =  0  ! Enter # of Secondary Outputs (not used).
      Ieqsol  =  1  ! Model type: 1 - DIFFEQ, 2 - AMAT, 3 - OUTPUT only.
      Descr   = ' REVER1  '

CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each System Parameter (eg. Psym(1)='Kel')         C
C----c-----------------------------------------------------------------C
		Psym(1)='CL12'
		Psym(2)='CL21'
		Psym(3)='CL10'
		Psym(4)='CL20'
		Psym(5)='V1'
		Psym(6)='V2'
		Psym(7)='V3'
		Psym(8)='CLD1'
		Psym(9)='CLD2'
		Psym(10)='V4'

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
		Real*8 CL12,CL21,CL10,CL20,V1,V2,CLD1,V3,CLD2,V4
C----------------------------------------------------------------------C
C   Enter Differential Equations Below  {e.g.  XP(1) = -P(1)*X(1) }    C
C----c-----------------------------------------------------------------C
		CL12=P(1)
		CL21=P(2)
		CL10=P(3)
		CL20=P(4)
		V1=P(5)
		V2=P(6)
		V3=P(7)
		CLD1=P(8)
		CLD2=P(9)
		V4=P(10)
		
		XP(1)=-(CL10+CL12+CLD1)*X(1)/V1+CL21*X(2)/V1+CLD1*X(5)/V1
		XP(2)=CL12*X(1)/V2-(CL20+CL21)*X(2)/V2+CLD2*X(6)/V2-CLD2*X(2)/V2
		XP(3)=-(CL20+CL21+CLD2)*X(3)/V2+CL12*X(4)/V2+CLD2*X(6)/V2
		XP(4)=CL21*X(3)/V1-(CL10+CL12)*X(4)/V1+CLD1*X(5)/V1-CLD1*X(4)/V1
		XP(5) = CLD1*X(1)/V3-CLD1*X(5)/V3
		XP(6) = CLD2*X(3)/V4-CLD2*X(6)/V4
		XP(7) = CLD1*X(4)/V3-CLD1*X(5)/V3
		XP(8) = CLD2*X(2)/V4-CLD1*X(6)/V4
		
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

        Real*8 Y(MaxNOE),T,X(MaxNDE)

CC
C----------------------------------------------------------------------C
C   Enter Output Equations Below   {e.g.  Y(1) = X(1)/P(2) }           C
C----c-----------------------------------------------------------------C

        Y(1) = X(1)
		Y(2) = X(2)
		Y(3) = X(3)
		Y(4) = X(4)

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
		 V(3) = (PV(1) + PV(2)*Y(1))**2 
		 V(4) = (PV(1) + PV(2)*Y(2))**2

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