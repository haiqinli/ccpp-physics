!>\file gsd_chem_emis_driver.F90
!! This file is CCPP GSD Chemistry emission driver to provide water/ice friendely aerosols for Thmopson MP 
!! Haiqin.Li@noaa.gov 11/2020

 module gsd_chem_emis_driver

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use gsd_chem_config
   use dust_data_mod
   use seas_mod,        only : gocart_seasalt_driver
   use dust_fengsha_mod,only : gocart_dust_fengsha_driver
   use opt_mod

   implicit none

   private

   public :: gsd_chem_emis_driver_init, gsd_chem_emis_driver_run, gsd_chem_emis_driver_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine gsd_chem_emis_driver_init()
      end subroutine gsd_chem_emis_driver_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_gsd_chem_emis_driver_finalize Argument Table
!!
      subroutine gsd_chem_emis_driver_finalize()
      end subroutine gsd_chem_emis_driver_finalize

!> \defgroup gsd_chem_group GSD Chem emission driver Module
!! This is the gsd chemistry
!>\defgroup gsd_chem_emis_driver GSD Chem emission driver Module  
!> \ingroup gsd_chem_group
!! This is the GSD Chem emission driver Module
!! \section arg_table_gsd_chem_emis_driver_run Argument Table
!! \htmlinclude gsd_chem_emis_driver_run.html
!!
!>\section gsd_chem_emis_driver GSD Chemistry Scheme General Algorithm
!> @{
    subroutine gsd_chem_emis_driver_run(im, kte, kme, ktau, dt, garea, land,           &
                   u10m, v10m, ustar, rlat, rlon, tskin, pb2d,                         &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum, w,              &
                   nsoil, smc, vegtype, soiltyp, sigmaf, dswsfc, zorl,snow,            &
                   dust_in, emi_in, fire_GBBEPx, ntrac, gq0,tile_num,                  &
                   imp_physics, imp_physics_thompson, nwfa, nifa, nwem2d, niem2d,      &
                   emdust, emseas, emsulf, emanoc, emfroc, maod,                       &
                   errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau,nsoil,tile_num
    integer,        intent(in) :: ntrac
    real(kind_phys),intent(in) :: dt

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    integer, dimension(im), intent(in) :: land, vegtype, soiltyp        
    real(kind_phys), dimension(im,nsoil), intent(in) :: smc
    real(kind_phys), dimension(im,    5), intent(in) :: dust_in
    real(kind_phys), dimension(im,   10), intent(in) :: emi_in
    real(kind_phys), dimension(im,    5), intent(in) :: fire_GBBEPx
    real(kind_phys), dimension(im), intent(in) :: u10m, v10m, ustar,              &
                garea, rlat,rlon, tskin, pb2d, sigmaf, dswsfc, zorl, snow
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d,         &
                us3d, vs3d, spechum, w
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0
    real(kind_phys), dimension(im    ), intent(inout) :: maod
    real(kind_phys), dimension(im,kte), intent(inout) :: nwfa, nifa
    real(kind_phys), dimension(im    ), intent(inout) :: nwem2d, niem2d
    real(kind_phys), dimension(im    ), intent(inout) :: emdust, emseas, emsulf, emanoc, emfroc
    integer, intent(in   ) :: imp_physics, imp_physics_thompson
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy, u_phy, v_phy,  &
                     p_phy, z_at_w, dz8w, p8w, t8w, rho_phy, vvel, zmid

    real(kind_phys), dimension(ims:im, jms:jme) :: u10, v10, ust, tsk,            &
                     xland, xlat, xlong, dxy, pbl

!>- sea salt & chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    real(kind_phys), dimension(ims:im, 1, jms:jme, 1:num_emis_seas  )  :: emis_seas
    real(kind_phys), dimension(ims:im, jms:jme) :: seashelp

    integer :: ide, ime, ite, kde

!>- dust & chemistry variables
    real(kind_phys), dimension(ims:im, jms:jme) :: ssm, rdrag, uthr, snowh  ! fengsha dust
    real(kind_phys), dimension(ims:im, jms:jme) :: vegfrac, rmol, gsw, znt, clayf, sandf
    real(kind_phys), dimension(ims:im, nsoil, jms:jme) :: smois
    real(kind_phys), dimension(ims:im, 1:1, jms:jme, 1:num_emis_dust) :: emis_dust
    integer,         dimension(ims:im, jms:jme) :: isltyp, ivgtyp

!>- optical variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: relhum
    real(kind_phys), dimension(ims:im,          jms:jme) :: aod
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:nbands) :: extt, ssca, asympar

    ! -- buffers
    real(kind_phys), dimension(ims:im, jms:jme, num_ebu_in) :: ebu_in
    real(kind_phys), dimension(ims:im, kms:kemit, jms:jme, 1:num_emis_ant) :: emis_ant

    ! -- parameter to caluclate wfa&ifa (m)
    !real(kind_phys), parameter :: mean_diameter1= 4.E-8, sigma1=1.8
    !real(kind_phys), parameter :: mean_diameter2= 4.E-7, sigma2=1.8
    !!real(kind_phys), parameter :: mean_diameter1= 5.E-8, sigma1=1.8
    real(kind_phys), parameter :: mean_diameter1= 6.E-8, sigma1=1.8
    real(kind_phys), parameter :: mean_diameter2= 12.E-7, sigma2=1.8
    !real(kind_phys), parameter :: mean_diameter2= 6.E-7, sigma2=1.8
    !-- aerosol density (kg/m3)
    real(kind_phys), parameter :: density_dust= 2.6e+3, density_sulfate=1.8e+3
    real(kind_phys), parameter :: density_oc  = 1.4e+3, density_seasalt=2.2e+3
    real(kind_phys), dimension(ndust), parameter :: den_dust  = (/2500., 2650., 2650., 2650., 2650./)
    real(kind_phys), dimension(im) :: emis_sulfate, emis_oc, emis_seasalt, emis_dusts
    real(kind_phys), dimension(im) :: emis_anoc, emis_fireoc
    real(kind_phys), dimension(im) :: daero_emis_wfa, daero_emis_ifa
    real(kind_phys) :: fact_wfa, fact_ifa


!>-- local variables
    logical :: call_radiation
    logical :: store_arrays
    integer :: nbegin, nv, nvv
    integer :: i, j, jp, k, kp, n
  

    errmsg = ''
    errflg = 0

!    print*,'hli ktau',ktau
    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    emis_sulfate   = 0.
    emis_oc        = 0.
    emis_seas      = 0.
    emis_seasalt   = 0.
    emis_dust      = 0.
    emis_dusts     = 0.
    emis_anoc      = 0.
    emis_fireoc    = 0.
    fact_wfa       = 0.
    fact_ifa       = 0.
    daero_emis_wfa = 0.
    daero_emis_ifa = 0.



!>- get ready for chemistry run
    call chem_emis_prep(                                            &
        ktau,                      &
        u10m,v10m,ustar,land,garea,rlat,rlon,tskin,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                 &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                   &
        snow,dust_in,emi_in,fire_GBBEPx,                            &
        pb2d,                               &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        t8w,                                               &
        z_at_w,vvel,zmid,                                               &
        ntrac,gq0,                                                      &
        num_chem, num_moist, num_ebu_in,                                &
        num_emis_ant,                     &
        emis_ant,                                             &
        moist,chem,relhum,ebu_in,                                    &
        smois,ivgtyp,isltyp,vegfrac,rmol,gsw,znt,pbl,               &
        snowh,clayf,rdrag,sandf,ssm,uthr,             &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)

!>- compute sea-salt
    ! -- compute sea salt
    if (seas_opt >= SEAS_OPT_DEFAULT) then
    call gocart_seasalt_driver(ktau,dt,rri,t_phy,moist,                 &
        u_phy,v_phy,chem,rho_phy,dz8w,u10,v10,ust,p8w,tsk,              &
        xland,xlat,xlong,dxy,g,emis_seas,                               &
        seashelp,num_emis_seas,num_moist,num_chem,seas_opt,             &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)
    endif
    !print*,'hli emis_seas',maxval(emis_seas(:,1,1,3))

    !-- compute dust
    select case (dust_opt)
      case (DUST_OPT_FENGSHA)
       dust_alpha     = 2.0
       dust_gamma     = 1.8
       call gocart_dust_fengsha_driver(dt,chem,rho_phy,smois,p8w,ssm,   &
            isltyp,vegfrac,snowh,xland,dxy,g,emis_dust,ust,znt,         &
            clayf,sandf,rdrag,uthr,                                     &
            num_emis_dust,num_moist,num_chem,nsoil,                     &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte)
    end select
    !print*,'hli emis_dust',maxval(emis_dust(:,1,1,1))


    

   if (imp_physics == imp_physics_thompson) then

    do i = 1, im
     !do n = 1, num_emis_dust
     ! emis_dusts(i) = emis_dusts(i)+ emis_dust(i,1,1,n) ! dust emission :ug/m2/s
     !enddo
     !do n = 1, num_emis_seas
     ! emis_seasalt(i) = emis_seasalt(i) + emis_seas(i,1,1,n)*1.e+9 ! sea salt emission: ug/m2/s 
     !enddo
      emis_dusts  (i) = emis_dusts  (i) + emis_dust(i,1,1,1) ! dust emission :ug/m2/s (bin 1 dust)
      emis_seasalt(i) = emis_seasalt(i) + emis_seas(i,1,1,1)*1.e+9 ! sea salt emission: ug/m2/s  (bin 1 sea salt)
    
     emis_sulfate(i) = emis_ant(i,kts,1,p_e_sulf) ! anthro. sulfate, normally it is zero in HTAP and CEDS
     emis_anoc   (i) = emis_ant(i,kts,1,p_e_oc  ) ! anthro. organic carbon: ug/(m2*s)
     emis_fireoc (i) = ebu_in  (i,kts,p_e_oc  )   ! fire organic carbon: ug/(m2*s)
     emis_oc     (i) = (emis_anoc(i) + emis_fireoc(i))*0.125
     !emis_oc     (i) = emis_anoc(i) + emis_fireoc(i)
     ! hli output emission for diagnostic purpose
     emdust(i)=emis_dusts  (i)
     emseas(i)=emis_seasalt(i)
     emsulf(i)=emis_sulfate(i)
     emanoc(i)=emis_anoc   (i)
     emfroc(i)=emis_fireoc (i)
    enddo
    !print*,'hli ktau, emis_fireoc',ktau,maxval(emis_fireoc)

    fact_wfa = 1.e-9*6.0/pi*exp(4.5*log(sigma1)**2)/mean_diameter1**3 
    fact_ifa = 1.e-9*6.0/pi*exp(4.5*log(sigma2)**2)/mean_diameter2**3 
    do i = 1, im
     daero_emis_wfa(i) = emis_sulfate(i)/density_sulfate + emis_oc(i)/density_oc   &
                       + emis_seasalt(i)/density_seasalt
     daero_emis_ifa(i) = emis_dusts(i)/density_dust

     daero_emis_wfa(i) = daero_emis_wfa(i)*fact_wfa*rri(i,kemit,1)/dz8w(i,kemit,1)
     daero_emis_ifa(i) = daero_emis_ifa(i)*fact_ifa*rri(i,kemit,1)/dz8w(i,kemit,1)

     nwem2d(i)=daero_emis_wfa(i)
     niem2d(i)=daero_emis_ifa(i)

     nwfa(i,kemit)     = nwfa(i,kemit) + daero_emis_wfa(i)*dt
     nifa(i,kemit)     = nifa(i,kemit) + daero_emis_ifa(i)*dt
    enddo
     !print*,'hli rri',maxval(rri(:,kemit,1))
     !print*,'hli dz8w',maxval(dz8w(:,kemit,1))
     !print*,'hli daero_emis_wfa,nwfa',maxval(daero_emis_wfa),maxval(nwfa(:,kemit))
     !if(maxval(daero_emis_ifa).gt.0.)then
     !print*,'hli daero_emis_ifa,nifa',maxval(daero_emis_ifa),maxval(nifa(:,kemit))
     !endif


!>-- compuate AOD ---
    aer_ra_feedback = 2
    if (call_radiation) then
      store_arrays = .false.
      select case (aer_ra_feedback)
        case (2)
!          call aero_opt('sw',dz8w,chem,nwfa,nifa                        &
!                   ,rri,relhum,aod                                      &
!!                  ,extt,ssca,asympar                                   &
!                   ,extt,ssca,asympar,num_chem                          &
!                   ,ids,ide, jds,jde, kds,kde                           &
!                   ,ims,ime, jms,jme, kms,kme                           &
!                   ,its,ite, jts,jte, kts,kte)
          store_arrays = .true.
        case default
          ! -- no feedback
      end select
      if (store_arrays) then
        maod(its:ite) = aod(its:ite,1)
      end if
    endif

   endif
!
 end subroutine gsd_chem_emis_driver_run

 subroutine chem_emis_prep(                                        &
        ktau,                     &
        u10m,v10m,ustar,land,garea,rlat,rlon,ts2d,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                  &
        snow_cpl,dust_in,emi_in,fire_GBBEPx,                           &
        pb2d,                              &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                          &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                  &
        t8w,                                              &
        z_at_w,vvel,zmid,                                              &
        ntrac,gq0,                                                     &
        num_chem, num_moist, num_ebu_in,                               &
        num_emis_ant,                    &
        emis_ant,                                             &
        moist,chem,relhum,ebu_in,                                   &
        smois,ivgtyp,isltyp,vegfrac,rmol,gsw,znt,pbl,              &
        snowh,clayf,rdrag,sandf,ssm,uthr,            &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    integer, intent(in) :: ktau

    !FV3 input variables
    integer, intent(in) :: nsoil
    integer, dimension(ims:ime), intent(in) :: land, vegtype, soiltyp
    integer, intent(in) :: ntrac
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         u10m, v10m, ustar, garea, rlat, rlon, ts2d, sigmaf, dswsfc,       &
         zorl, snow_cpl, pb2d
    real(kind=kind_phys), dimension(ims:ime, nsoil),   intent(in) :: smc 
    real(kind=kind_phys), dimension(ims:ime,     5),   intent(in) :: dust_in
    real(kind=kind_phys), dimension(ims:ime,    10),   intent(in) :: emi_in
    real(kind=kind_phys), dimension(ims:ime,     5),   intent(in) :: fire_GBBEPx
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,w
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_moist, num_emis_ant, num_ebu_in
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte


    real(kind_phys), dimension(ims:ime, kms:kemit, jms:jme, num_emis_ant), intent(inout) :: emis_ant
    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in),intent(out) :: ebu_in
    
    integer,dimension(ims:ime, jms:jme), intent(out) :: isltyp, ivgtyp
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, relhum, t_phy, u_phy, v_phy, p_phy, rho_phy, dz8w, p8w, t8w, vvel, zmid
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         u10, v10, ust, tsk, xland, xlat, xlong, dxy, vegfrac, rmol, gsw, znt,     &
         pbl, snowh, clayf, rdrag, sandf, ssm, uthr
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w
    real(kind_phys), dimension(ims:ime, nsoil, jms:jme), intent(out) :: smois
    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in) :: emiss_ab
    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in) :: emiss_abu
    real(kind_phys), parameter :: frac_so2_ant = 0.5_kind_phys  ! antropogenic so2 fraction
    real(kind_phys), parameter :: frpc  = 1.e+09_kind_phys      ! FRP conversion factor

    ! -- local variables
    integer i,ip,j,jp,k,kp,kk,kkp,nv,l,ll,n


    ! -- initialize output arrays
    isltyp         = 0._kind_phys
    ivgtyp         = 0._kind_phys
    rri            = 0._kind_phys
    relhum         = 0._kind_phys
    t_phy          = 0._kind_phys
    u_phy          = 0._kind_phys
    v_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    t8w            = 0._kind_phys
    vvel           = 0._kind_phys
    zmid           = 0._kind_phys
    u10            = 0._kind_phys
    v10            = 0._kind_phys
    ust            = 0._kind_phys
    tsk            = 0._kind_phys
    xland          = 0._kind_phys
    xlat           = 0._kind_phys
    xlong          = 0._kind_phys
    dxy            = 0._kind_phys
    vegfrac        = 0._kind_phys
    rmol           = 0._kind_phys
    gsw            = 0._kind_phys
    znt            = 0._kind_phys
    pbl            = 0._kind_phys
    snowh          = 0._kind_phys
    clayf          = 0._kind_phys
    rdrag          = 0._kind_phys
    sandf          = 0._kind_phys 
    ssm            = 0._kind_phys
    uthr           = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys

    chem = 0._kind_phys

    do i=its,ite
     u10  (i,1)=u10m (i)
     v10  (i,1)=v10m (i)
     tsk  (i,1)=ts2d (i)
     ust  (i,1)=ustar(i)
     dxy  (i,1)=garea(i)
     xland(i,1)=real(land(i))
     xlat (i,1)=rlat(i)*180./pi
     xlong(i,1)=rlon(i)*180./pi
     gsw  (i,1)=dswsfc(i)
     znt  (i,1)=zorl(i)*0.01
     pbl  (i,1)=pb2d(i)
     snowh(i,1)=snow_cpl(i)*0.001
     clayf(i,1)=dust_in(i,1)
     rdrag(i,1)=dust_in(i,2)
     sandf(i,1)=dust_in(i,3)
     ssm  (i,1)=dust_in(i,4)
     uthr (i,1)=dust_in(i,5)
     ivgtyp (i,1)=vegtype(i)
     isltyp (i,1)=soiltyp(i)
     vegfrac(i,1)=sigmaf (i)
    enddo
   
    rmol=0.

    do k=1,nsoil
     do j=jts,jte
      do i=its,ite
       smois(i,k,j)=smc(i,k)
      enddo
     enddo
    enddo

    if (ktau <= 1) then
      emis_ant = 0.
     !emis_vol = 0.
    end if

    do j=jts,jte
      jp = j - jts + 1
      do i=its,ite
         ip = i - its + 1
         z_at_w(i,kts,j)=max(0.,ph3d(ip,1)/g)
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=abs(ph3d(ip,kp+1)-ph3d(ip,kp))/g
          z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          p8w(i,k,j)=pr3d(ip,kp)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kk=min(k,kte)
        kkp = kk - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
          t_phy(i,k,j)=tk3d(ip,kkp)
          p_phy(i,k,j)=prl3d(ip,kkp)
          u_phy(i,k,j)=us3d(ip,kkp)
          v_phy(i,k,j)=vs3d(ip,kkp)
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          vvel(i,k,j)=-w(ip,kkp)*rri(i,k,j)/g 
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=gq0(ip,kkp,p_atm_shum)
          if (t_phy(i,k,j) > 265.) then
            moist(i,k,j,2)=gq0(ip,kkp,p_atm_cldq)
            moist(i,k,j,3)=0.
            if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          else
            moist(i,k,j,2)=0.
            moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldq)
            if(moist(i,k,j,3) < 1.e-8)moist(i,k,j,3)=0.
          endif
          relhum(i,k,j) = .95
          relhum(i,k,j) = MIN( .95, moist(i,k,j,1) / &
            (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
            (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
          relhum(i,k,j)=max(0.1,relhum(i,k,j))
          !--
          zmid(i,k,j)=phl3d(ip,kkp)/g
        enddo
      enddo
    enddo

    do j=jts,jte
      do k=2,kte
        do i=its,ite
          t8w(i,k,j)=.5*(t_phy(i,k,j)+t_phy(i,k-1,j))
        enddo
      enddo
    enddo

    ! -- only used in phtolysis....
    do j=jts,jte
      do i=its,ite
        t8w(i,1,j)=t_phy(i,1,j)
        t8w(i,kte+1,j)=t_phy(i,kte,j)
      enddo
    enddo

    ! -- fire
    emiss_ab  = 0.   ! background
    emiss_abu = 0.   ! background
    do j=jts,jte
     do i=its,ite
      emiss_ab(i,j,p_e_bc)   =emi_in(i,1)
      emiss_ab(i,j,p_e_oc)   =emi_in(i,2)
      emiss_ab(i,j,p_e_sulf) =emi_in(i,3)
      emiss_ab(i,j,p_e_pm_25)=emi_in(i,4)
      emiss_ab(i,j,p_e_so2)  =emi_in(i,5)
      emiss_ab(i,j,p_e_pm_10)=emi_in(i,6)
     enddo
    enddo

    select case (plumerise_flag)
      case (FIRE_OPT_GBBEPx)
        do j=jts,jte
         do i=its,ite
          emiss_abu(i,j,p_e_bc)   =fire_GBBEPx(i,1)
          emiss_abu(i,j,p_e_oc)   =fire_GBBEPx(i,2)
          emiss_abu(i,j,p_e_pm_25)=fire_GBBEPx(i,3)
          emiss_abu(i,j,p_e_so2)  =fire_GBBEPx(i,4)
         !plume(i,j,1)            =fire_GBBEPx(i,5)
         enddo
        enddo
      case default
          ! -- no further option available
    end select


    k=kts
    if (p_bc2 > 1) then
      do j=jts,jte
        do i=its,ite
          emis_ant(i,k,j,p_e_bc)=emiss_ab(i,j,p_e_bc)
          emis_ant(i,k,j,p_e_oc)=emiss_ab(i,j,p_e_oc) 
         !emis_ant(i,k,j,p_e_oc)=emiss_ab(i,j,p_e_oc) + emiss_ab(i,j,p_e_pm_25)
          emis_ant(i,k,j,p_e_sulf)=emiss_ab(i,j,p_e_sulf)
          emis_ant(i,k,j,p_e_so2)=frac_so2_ant * emiss_ab(i,j,p_e_so2)
          emis_ant(i,k,j,p_e_dms)= 0. !emiss_ab(j,p_e_dms)
          emis_ant(i,k,j,p_e_pm_25)=emiss_ab(i,j,p_e_pm_25)
          emis_ant(i,k,j,p_e_pm_10)=emiss_ab(i,j,p_e_pm_10)

          ebu_in(i,j,p_ebu_in_pm10)=emiss_abu(i,j,p_e_pm_10)
          ebu_in(i,j,p_ebu_in_dms)= 0._kind_phys

          select case (plumerise_flag)
            case (FIRE_OPT_GBBEPx)
              ebu_in(i,j,p_ebu_in_oc)   = frpc * (emiss_abu(i,j,p_e_oc))
             !ebu_in(i,j,p_ebu_in_oc)   = frpc * (emiss_abu(i,j,p_e_pm_25) - emiss_abu(i,j,p_e_bc))
              ebu_in(i,j,p_ebu_in_bc)   = frpc * emiss_abu(i,j,p_e_bc)
              ebu_in(i,j,p_ebu_in_pm25) = frpc * (emiss_abu(i,j,p_e_pm_25) - emiss_abu(i,j,p_e_bc) - emiss_abu(i,j,p_e_oc))
              ebu_in(i,j,p_ebu_in_so2)  = frpc * emiss_abu(i,j,p_e_so2)
            case default
              ! -- no further option available
          end select
        enddo
      enddo
    endif
 


  end subroutine chem_emis_prep

!> @}
  end module gsd_chem_emis_driver
