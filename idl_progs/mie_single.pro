PRO Mie_single,Dx,Cm,Dqv=dqv,Dqxt,Dqsc,Dqbk,Dg,Xs1,Xs2,Dph,dlm=dlm
;+
; NAME:
;       MIE_SINGLE
; PURPOSE:
;       Calculates the scattering parameters of a series of particles
;       using the Mie scattering theory.
; CATEGORY:
;       EODG Mie routines
; CALLING SEQUENCE:
;       mie_single, Dx, Cm, Inp [, Dqv = dqv] $
;       [, Dqxt] [, Dqsc] [, Dqbk] [, Dg] [, Xs1] [, Xs2] [, Dph]
; INPUTS:
;       Dx:         A 1D array of particle size parameters
;       Cm:         The complex refractive index of the particles
;       Inp:        Number of scattering angles at which to calculate
;                   intensity functions etc
; OPTIONAL KEYWORD PARAMETERS:
;       dqv:        An array of the cosines of scattering angles at
;                   which to compute the phase function.
;       dlm:        If set the code will call the IDL DLM version of the
;                   algorithm
; OUTPUT PARAMETERS:
;       Dqxt:       The extinction efficiency
;       Dqsc:       The scattering efficiency
;       Dqbk:       The backscattering efficiency
;       Dg:         The asymmetry parameter
;       Xs1:        The first amplitude function - amplitude of light
;                   polarized in the plane perpendicular to the
;                   directions of incident light propagation and
;                   observation.
;       Xs2:        The second amplitude function - amplitude of light
;                   polarized in the plane parallel to the directions
;                   of incident light propagation and observation.
;                   NB. Xs1 and Xs2 are complex arrays of the same
;                   dimension as Dqv and are only calculated if Dqv is
;                   specified.
;       Dph:        The phase function - an array of the same
;                   dimension as Dqv. Also only calculated if Dqv is
;                   specified.
; MODIFICATION HISTORY
;       G. Thomas 1998 mie_uoc.pro (translation of mieint.f to IDL)
;
;       D. Grainger 2001(?) mie_uoc_d.pro (Added support for arrays of
;       particle sizes and included calculation of phase function)
;
;       G. Thomas Sept. 2003 (Put into EODG routines format)
;
;       G. Thomas Feb. 2004 (Introduced explicit double precision
;       numerical values into all computational expressions)
;       G. Thomas Apr. 2005 (NmX assignment changed to type long)
;       G. Thomas Apr. 2005 (Added dlm keyword)
;       G. Thomas Apr. 2005 (Changed code to ensure Qbsc is always
;       caculated for the backscatter direction)
;       G. Thomas Jun. 2005 (Added calculation of phase function after
;       calling mie_dlm_single, since the DLM no longer returns it)
;       G. Thomas Nov. 2006 (Changed the calculation of the 
;       backscatter efficiency to be done directly from the A and B
;       values, rather than from the intensity at 180 degrees. Also, 
;       fixed a small bug with passing Dqv to the DLM)
;       A. Smith 12 Oct 07  (Fixed bug in non-DLM calculation of Dph)
;       G. Thomas Jun. 2007 (Bug fix in calculation of backscatter
;                            efficiency)
;-
;   If an array of cos(theta) is provided, calculate phase function

  Imaxx = 12000l
  RIMax = 2.5
  Itermax = long(Imaxx * RIMax)
  Imaxnp = 10000l     ; Change this as required
  Sizes = N_Elements(Dx)

  IF KEYWORD_SET(dlm) THEN BEGIN
  ; If the dlm keyword is set, use the DLM version of the code
  DxArr = dblarr(Sizes)
  DxArr[*] = Dx
  DCm = dcomplex(Cm)
  IF KEYWORD_SET(dqv) THEN BEGIN
      Inp = n_elements(dqv)
      mie_dlm_single,DxArr,DCm,dqv=double(dqv),Dqxt,Dqsc,Dqbk,Dg,Xs1,Xs2
      ; Cannot get the DLM to return the phase function. Very mysterious...
            Dph = dblarr(Inp,Sizes)
            for i = 0,Sizes-1 do $
               Dph[*,i] = 2d0 * double(Xs1[*,i]*CONJ(Xs1[*,i]) + $
                                       Xs2[*,i]*CONJ(Xs2[*,i])) $
                          / (DxArr[i]^2 * Dqsc[i])
  ENDIF ELSE BEGIN
      mie_dlm_single,DxArr,DCm,Dqxt,Dqsc,Dqbk,Dg
  ENDELSE

  ENDIF ELSE BEGIN ;No DLM? Do everything in IDL

  if KEYWORD_SET(dqv) then begin
     tmp = where(dqv eq -1.0,bktheta)
     if bktheta eq 0 then dqv2 = [dqv,-1d] else dqv2 = dqv
     Inp  = n_elements(dqv)
     Inp2 = n_elements(dqv2)
     ph = dblarr(Inp)
     Xs1 = ComplexArr(Inp,Sizes)
     Xs2 = ComplexArr(Inp,Sizes)
     Dph = Dblarr(Inp,Sizes)
  endif else begin
      Inp = 1
      Inp2 = Inp
      Dqv = [-1d0]
      Dqv2 = Dqv
  endelse

  Dqxt = Dblarr(Sizes)
  Dqsc = Dblarr(Sizes)
  Dqbk = DComplexArr(Sizes)
  Dg = Dblarr(Sizes)

  For Size = 0l, Sizes - 1 Do Begin

    IF (Dx(Size) GT Imaxx) THEN MESSAGE, 'Error: Size Parameter Overflow in Mie'
    Ir = 1.D0 / Cm
    Y =  Dx(Size) * Cm

    IF (Dx(Size) LT 0.02) THEN  NStop = 2 ELSE BEGIN
      CASE 1 OF
        (Dx(Size) LE 8.0)    : NStop = Dx(Size) + 4.00*Dx(Size)^(1./3.) + 2.0
        (Dx(Size) LT 4200.0) : NStop = Dx(Size) + 4.05*Dx(Size)^(1./3.) + 2.0
        ELSE                 : NStop = Dx(Size) + 4.00*Dx(Size)^(1./3.) + 2.0
      ENDCASE
    END
    NmX = LONG(MAX([NStop,ABS(Y)]) + 15.)
    D = DCOMPLEXARR(Nmx+1)

    FOR N = Nmx-1,1,-1 DO BEGIN
      A1 = (N+1) / Y
      D(N) = A1 - 1/(A1+D(N+1))
    END

    Sm = DCOMPLEXARR(Inp2)
    Sp = DCOMPLEXARR(Inp2)
    Pi0 = DCOMPLEXARR(Inp2)
    Pi1 = DCOMPLEX(REPLICATE(1.D0,Inp2),REPLICATE(0.D0,Inp2))

    Psi0 = Cos(Dx(Size))
    Psi1 = Sin(Dx(Size))
    Chi0 =-Sin(Dx(Size))
    Chi1 = Cos(Dx(Size))
    Xi0 = DCOMPLEX(Psi0,Chi0)
    Xi1 = DCOMPLEX(Psi1,Chi1)

    Dg(Size) = 0.D0
    Dqsc(Size) = 0.D0
    Dqxt(Size) = 0.D0
    Dqbk(Size) = 0.D0
    Tnp1 = 1D0

    FOR N = 1l,Nstop DO BEGIN
      DN = Double(N)
      Tnp1 = Tnp1 + 2D0
      Tnm1 = Tnp1 - 2D0
      A2 = Tnp1 / (DN*(DN+1.D0))
      Turbo = (DN+1.D0) / DN
      Rnx = DN/Dx(Size)
      Psi = DOUBLE(Tnm1)*Psi1/Dx(Size) - Psi0
      Chi = Tnm1*Chi1/Dx(Size)       - Chi0
      Xi = DCOMPLEX(Psi,Chi)
      A = ((D[N]*Ir+Rnx)*Psi-Psi1) / ((D[N]*Ir+Rnx)*  Xi-  Xi1)
      B = ((D[N]*Cm+Rnx)*Psi-Psi1) / ((D[N]*Cm+Rnx)*  Xi-  Xi1)
      Dqxt(Size) = Tnp1 * DOUBLE(A + B)                 + Dqxt(Size)
      Dqsc(Size) = Tnp1 * DOUBLE(A*CONJ(A) + B*CONJ(B)) + Dqsc(Size)
      IF (N GT 1) THEN Dg(Size) = Dg(Size) $
                + (dN*dN - 1) * DOUBLE(ANM1*CONJ(A) + BNM1 * CONJ(B)) / dN $
                + TNM1 * DOUBLE(ANM1*CONJ(BNM1)) / (dN*dN - dN)
      Anm1 = A
      Bnm1 = B


      S = Dqv2 * Pi1
      T = S - Pi0
      if arg_present(Dqbk) or n_elements(dph) gt 0 then begin
         Taun = N*T - Pi0
         Sp = (A2 * (A + B)) * (Pi1 + Taun) + Sp
         Sm = (A2 * (A - B)) * (Pi1 - Taun) + Sm
      endif
      Pi0 = Pi1
      Pi1 = S + T*Turbo

      Psi0 = Psi1
      Psi1 = Psi
      Chi0 = Chi1
      Chi1 = Chi
      Xi1 = DCOMPLEX(Psi1,Chi1)

    END; For Nstop

    IF (Dg(Size) GT 0) THEN Dg(Size) = 2D0 * Dg(Size) / Dqsc(Size)
;   The following lines are not needed unless dqv was set
    IF n_elements(dph) gt 0 THEN BEGIN  
        Xs1(*,Size) = ((Sp[0:Inp-1] + Sm[0:Inp-1]) / 2D0)
        Xs2(*,Size) = ((Sp[0:Inp-1] - Sm[0:Inp-1]) / 2D0)
        Dph(*,Size) = DOUBLE(Xs1(*,Size)*CONJ(Xs1(*,Size)) + $
                             Xs2(*,Size)*CONJ(Xs2(*,Size))) $
                            / Dqsc(Size)
    ENDIF

    if arg_present(Dqbk) then $
       Dqbk(Size) = ( Sp[Inp2-1] + Sm[Inp2-1] ) / 2d0     + Dqbk(Size)

  EndFor ; END of size loop

  Dqsc =  2D0 * Dqsc / Dx^2
  Dqxt =  2D0 * Dqxt / Dx^2
  if arg_present(Dqbk) then $
     Dqbk =  DOUBLE(Dqbk*CONJ(Dqbk)) / Dx^2 / !dpi


  ENDELSE ; END of if DLM keyword

End
