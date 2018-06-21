      program visr
c --------------------------------------------------------------------
c   visr (velocity interpolation to strain rates): a program to interpolate horizontal station velocities 
c   into on spot translation, rotation, and strain rates.  
c
c   For algorithm of the interpolation, see Shen et al., Optimal Interpolation of Spatially Discretized 
c   Geodetic Data, Bull. Seismol. Soc. Am., 2015. 
c
c   to run the program:
c   % visr < drive.file
c
c ---------------------------------------------------------------------
      implicit real*8 (a-h,l,o-z)
      parameter (maxsit=5000,maxc=10)
      character*100 head
      character*20 stnf,sitf,outf,crpf
      character*8 stnl(maxsit+1)
      integer*8 nstn,nslct,ncrp,nsite,m
      dimension lonl(maxsit+1),latl(maxsit+1)
      dimension lon(maxsit),lat(maxsit),uxl(maxsit),sxl(maxsit),
     .          uyl(maxsit),syl(maxsit),cxy(maxsit)
      dimension pos(2,maxsit),area(maxsit),carea(maxsit)
      dimension alon(maxc),alat(maxc),blon(maxc),blat(maxc),
     .          dcs(maxc),dsn(maxc)
      dimension bb(6),ainv(6,6),v(6),x(3),covx(3,3),tx(1),ty(1)
      common /veldat/nstn,lon,lat,uxl,sxl,uyl,syl,cxy,area
      common /crpdat/ncrp,alon,alat,blon,blat,dcs,dsn
      common /inv_v/bb,ainv,rtau,wtt,nslct,chisq
      common /prmtr/ pi,cov,cutoff_dis,wt_az
*
*   setup the constants                          
*
      pi = 4.0d0*datan(1.0d0)
      cov = pi/180.0d0
*
*   read input file names and options 
*
      read(5,'(a)') stnf
      read(5,'(a)') outf
      read(5,*) is_wght
      write(*,111) 'is_wght', is_wght
      read(5,*) id_wght
      write(*,111) 'id_wght', id_wght
      read(5,*) min_tau, max_tau, in_tau
      write(*,112) 'min_d', min_tau, 'max_d', max_tau, 'd_step', in_tau
      read(5,*) wt0
      write(*,110) 'wt0', wt0
      read(5,*) rsga
      write(*,110) 'rsga', rsga
      read(5,*) id
      write(*,111) 'id', id
112   format(A,I7,A,I7,A,I5)
111   format(A,I5)
110   format(A,F7.2)

      if ( id.eq.1 ) then
        read(5,*) lonmin,lonmax,latmin,latmax
      endif

      if ( id.eq.2 ) then
        read(5,'(a)') sitf
        open(4,file=sitf,status='old')
        i = 1
10      read(4,20,end=15) stnl(1),lonl(1),latl(1)
        if ( i.eq.1 ) then
          latmin=latl(1)
          latmax=latl(1)
          lonmin=lonl(1)
          lonmax=lonl(1)
        else
          if (latl(1).lt.latmin) latmin=latl(1)
          if (latl(1).gt.latmax) latmax=latl(1)
          if (lonl(1).lt.lonmin) lonmin=lonl(1)
          if (lonl(1).gt.lonmax) lonmax=lonl(1)
        endif
        i = i + 1
        goto 10
15      close(4)
        nsite = i - 1
      endif
20    format(a8,2f10.4)

      if ( id.eq.3 ) then
        read(5,*) lonmin,lonmax,latmin,latmax,dlon,dlat
        write(*,113) 'region lon/lat',
     .  lonmin,lonmax,latmin,latmax,dlon,dlat
      endif
113   format(a15,6f7.2)

      read(5,*) ncrp
      if (ncrp.gt.maxc) then
        write(6,61) maxc
        stop
      endif
      if (ncrp .ne. 0) then
        read(5,'(a)') crpf
        open(4,file=crpf,status='old')
        do i = 1, ncrp
          read(4,*) alon(i),alat(i),blon(i),blat(i)
          if (alon(i).gt.180.0d0 ) alon(i) = alon(i)-360.0d0
          if (blon(i).gt.180.0d0 ) blon(i) = blon(i)-360.0d0
        enddo
        close(4)
      endif
*
*   setup constants for data reweighting according to options                          
*
      if (is_wght .eq. 1) cutoff_dis = 2.15d0
      if (is_wght .eq. 2) cutoff_dis = 10.0d0
      wt_az = 0.25d0                                                         ! relative weight of mean azimuth to individual azimuth
      np_site = 6                                                            ! number of points to compute mean distance used as diameters of circlar weighting area 
      cfa1 = 2.0d0                                                           ! coefficient for circlular area weighing 
      cfa2 = 2.0d0                                                           ! coefficient for circular area weighting as substitute of voronoi area weighting
*
*   setup constants for data reweighting base on options                          
*
      exlat = cutoff_dis*max_tau/110.0d0
      exlon = cutoff_dis*max_tau/110.0d0/dcos(latmax*cov)

      slatmin = latmin - exlat
      slatmax = latmax + exlat
      slonmin = lonmin - exlon
      slonmax = lonmax + exlon

      if ( slatmin.lt.-90.0d0 ) slatmin = -90.0d0
      if ( slatmax.gt.90.0d0 ) slatmax = 90.0d0

      if (lonmin.gt.180.0d0) lonmin = lonmin - 360.0d0
      if (lonmax.gt.180.0d0) lonmax = lonmax - 360.0d0
      if (slonmin.gt.180.0d0) slonmin = slonmin - 360.0d0
      if (slonmax.gt.180.0d0) slonmax = slonmax - 360.0d0
*
*   read solution file of station coordinate and velocity
*
      open(4,file=stnf,status='old')

      i = 1
30    read(4,'(a)',end=40) head
      if (head(1:1).eq.'*') goto 30
      read(head,35) stnl(i),lonl(i),latl(i),
     .         uxl(i),sxl(i),uyl(i),syl(i),cxy(i)
*     write(*,114) stnl(i),lonl(i),latl(i),uxl(i),sxl(i),
*    .  uyl(i),syl(i),cxy(i)
      if (lonl(i).gt.180.0d0 ) lonl(i) = lonl(i)-360.0d0
      if (lonl(i).lt.slonmin.or.lonl(i).gt.slonmax.
     .    or.latl(i).lt.slatmin.or.latl(i).gt.slatmax) goto 30
      if (sxl(i).lt.rsga) then
          write(*,1130) 'Corrected sigma from', sxl(i), 'to', rsga
          sxl(i)=rsga
      endif
      if (syl(i).lt.rsga) then
          write(*,1130) 'Corrected sigma from', syl(i), 'to', rsga
          syl(i)=rsga
      endif
*     write(*,114) stnl(i),lonl(i),latl(i),uxl(i),sxl(i),
*    .  uyl(i),syl(i),cxy(i)
      i = i + 1
      goto 30
35    format(a8,2f10.4,2(f7.2,f5.2),f7.3)
40    close(4)
      nstn=i-1
114   format(a5,2f10.4,2(f7.2,f5.2),f7.3)
      write(*,115) 'Done reading stations; # of sta is', nstn
115   format(a20,I5)
1130  format(a20,f7.2,a5,f7.2)
1131  format(a20,I5)

      iden = 0
      do i = 1, nstn-1
         do j = i+1, nstn
            if (latl(i).eq.latl(j) .and. lonl(i).eq.lonl(j)) then
               print*,i,stnl(i),j,stnl(j)
               iden = 1
            endif
         enddo
      enddo
*     if (iden .eq. 1) stop
*     write(*,*) 'No duplicates found in station list.'

      if (nstn.gt.maxsit) then
        write(6,60) maxsit
        stop
      endif

60    format('# of stations in input file > maxsit',i5,
     .   ', please increase maxsit ...')
61    format('# of creep faults > maxc',i3,
     .   ', please increase maxc ...')
*
*   compute mass center of network
*
      xmean=0.0d0
      ymean=0.0d0
      do i=1,nstn
        xmean=xmean+lonl(i)
        ymean=ymean+latl(i)
      enddo
      xmean=xmean/nstn
      ymean=ymean/nstn
*
*  convert geodetic coordinates into cartesian coordinates
*
      do i=1,nstn
        lon(i)=lonl(i)
        lat(i)=latl(i)
      enddo
* Τα lat και lon πλέον είναι οι τοποκεντρικές συν/νες των σταθμών, σε
* σύστημα με κέντρο το κέντρο-βάρους του δικτύου. Οι ελλιψοειδείς
* συν/νες παραμένουν στα lonl και latl σε μοίρες
      call llxy(ymean,xmean,lat,lon,nstn)
*
*  compute voronoi area 
*
      if ( id_wght.eq.2 ) then

        do i=1,nstn
          pos(1,i)=lon(i)
          pos(2,i)=lat(i)
        enddo

        call get_voronoi_area_version_init(nstn,pos,area)

        do i = 1,nstn
          call cmp_area1(nstn,pos,np_site,cfa1,i,carea(i))
          if ( area(i).eq.-1.0d0 ) area(i) = carea(i)
          if ( area(i).gt.cfa2*carea(i) ) area(i) = carea(i)
        enddo

      endif
*
*  convert geodetic coordinates into cartesian coordinates for creep fault and
*  compute creep fault's azimuth
*
      if ( ncrp.ne.0 ) then 
        call llxy(ymean,xmean,alat,alon,ncrp)
        call llxy(ymean,xmean,blat,blon,ncrp)
        do i = 1, ncrp
          dx = blon(i) - alon(i)
          dy = blat(i) - alat(i)
          ds = dsqrt(dx**2 + dy**2)
          dcs(i) = dx/ds
          dsn(i) = dy/ds
        enddo
      endif
*
*   start interpolation loop
*
      open(9,file=outf,status='unknown')
*
*   check velocity compatibility
*
      if ( id.eq.1 ) then 
        write(9,70) 
        write(9,71) 

        do i = 1,nstn

          ty(1) = lat(i)
          tx(1) = lon(i)
          call visr_core(ty(1),tx(1),min_tau,max_tau,in_tau,
     .                   is_wght,id_wght,wt0,id,indxx)

          if ( indxx.eq.0 ) then
            v(1)=dsqrt(ainv(1,1))
            v(2)=dsqrt(ainv(2,2))
            vcrrl=ainv(1,2)/v(1)/v(2)
            write(9,75) stnl(i),lonl(i),latl(i),
     .      uxl(i),sxl(i),uyl(i),syl(i),cxy(i),
     .      bb(1),v(1),bb(2),v(2),vcrrl,rtau,wtt,chisq,nslct
          else 
            write(9,76) stnl(i),lonl(i),latl(i),
     .      uxl(i),sxl(i),uyl(i),syl(i),cxy(i),indxx
          endif

        enddo
      endif

70    format('*station   longitude  latitude            observation',
     . 26x,'interpolation                dis.    weight    chisq',
     . 3x,'nsite')
71    format(35x,'vx     svx    vy     svy   cxy',
     . 7x,'vx     svx    vy     svy   cxy')
75    format(1x,a8,2f10.4,2(2x,4f7.1,f7.2),f9.1,f9.2,f9.2,i6)
76    format(1x,a8,2f10.4,2x,4f7.1,f7.2,2x,7hindex =,i2)
*
*   interpolate velocity
*
      if ( id.eq.2 ) then
        write(9,80)

        open(4,file=sitf,status='old')
 
        do j = 1, nsite
          read(4,81) stnl(nstn+1),lonl(nstn+1),latl(nstn+1)

          do i = 1,nstn
            if (lonl(nstn+1).eq.lonl(i).
     .      and.latl(nstn+1).eq.latl(i)) then
              write (9,85) stnl(nstn+1),lonl(nstn+1),latl(nstn+1),
     .        uxl(i),sxl(i),uyl(i),syl(i),cxy(i),0.0,0.0,0.0,0
              goto 82
            endif
          enddo
        
          ty(1) = latl(nstn+1)
          tx(1) = lonl(nstn+1)
          m = 1
          call llxy(ymean,xmean,ty,tx,m)
          call visr_core(ty(1),tx(1),min_tau,max_tau,in_tau,
     .                   is_wght,id_wght,wt0,id,indxx)

          if ( indxx.eq.0 ) then
            v(1)=dsqrt(ainv(1,1))
            v(2)=dsqrt(ainv(2,2))
            vcrrl=ainv(1,2)/v(1)/v(2)
            write(9,85) stnl(nstn+1),lonl(nstn+1),latl(nstn+1),
     .      bb(1),v(1),bb(2),v(2),vcrrl,rtau,wtt,chisq,nslct
          else
            write(9,86) stnl(nstn+1),lonl(nstn+1),latl(nstn+1),indxx
          endif

82      continue
        enddo 
        close(4)
      endif

80    format('*station   longitude  latitude     vx     svx    vy',
     . 5x,'svy   cxy      dis.    weight    chisq   nsite')
81    format(a8,2f10.4)
85    format(1x,a8,2f10.4,2x,4f7.1,f7.2,f9.1,2f9.2,i6)
86    format(1x,a8,2f10.4,2x,7hindex =, i2)
*
*   interpolate horizontal station velocities into on spot translation, rotation, and strain rates
*
      if ( id.eq.3 ) then
        write(9,90) 
        write(9,94) 
 
        dltlat = latmax - latmin
        dltlon = lonmax - lonmin
        nlat = dltlat/dlat + 1.01d0
        nlon = dltlon/dlon + 1.01d0

        do ilat = 1, nlat
          tlat = latmin + dlat*(ilat-1)
        do ilon = 1, nlon
          tlon = lonmin + dlon*(ilon-1)

          ty(1) = tlat
          tx(1) = tlon
          m = 1
          write(*,121) 'Computing strain at', tlon, tlat
121       format(a20,2f10.4)
          call llxy(ymean,xmean,ty,tx,m)
          call visr_core(ty(1),tx(1),min_tau,max_tau,in_tau,
     .                   is_wght,id_wght,wt0,id,indxx)

          write(*,2000) 'lon=', tlon, 'lat=', tlat
          write(*,2001) 'Ux= ', bb(1)
          write(*,2001) 'Uy= ', bb(2)
          write(*,2001) 'tx= ', bb(3)
          write(*,2001) 'txy=', bb(4)
          write(*,2001) 'ty= ', bb(5)
          write(*,2001) 'w=  ', bb(6)
2000      format(a5, f10.5, a5, f10.5)
2001      format(a5, f10.5)

          if ( indxx.eq.0) then
            do i = 1,6
              v(i) = dsqrt(ainv(i,i))
            enddo
            vcrrl=ainv(1,2)/v(1)/v(2)
          
            do i = 1,3
              x(i) = bb(2+i)
*             x(1)=tx, x(2)=txy, x(3)=ty
              do j = 1,3
                covx(i,j) = ainv(2+i,2+j)
              enddo
            enddo

            call cmp_strain(x,covx,emax,emaxsd,emin,eminsd,
     .           taumax,tausd,dexazim,azsd,dilat,ddilat)

            emax = emax*1000.0d0
            emaxsd = emaxsd*1000.0d0
            emin = emin*1000.0d0
            eminsd = eminsd*1000.0d0
            taumax = taumax*1000.0d0
            tausd = tausd*1000.0d0
            dilat = dilat*1000.0d0
            ddilat = ddilat*1000.0d0

            do i = 3,6
              bb(i) = bb(i)*1000.0d0
              v(i) = v(i)*1000.0d0
            enddo

            write(9,91) tlon,tlat,bb(1),v(1),bb(2),v(2),vcrrl,
     .      bb(6),v(6),(bb(j),v(j),j=3,5),
     .      emax,emaxsd,emin,eminsd,taumax,tausd,dexazim,azsd,
     .      dilat,ddilat,rtau,wtt,chisq,nslct
          else  
            write(9,92) tlon,tlat,indxx
          endif

        enddo
        enddo
      endif

90    format('longitude latitude      vx+dvx       vy+dvy    cxy',
     . 9x,'w+dw         exx+dexx         exy+dexy         eyy+deyy',
     . 9x,'emax+demax       emin+demin        shr+dshr         azi+dazi',
     . 8x,'dilat+ddilat     dis.    weight    chisq   nsite')
94    format('  deg       deg          mm/a         mm/a      /'
     . 6x,'nano-rad/a    nano-strain/a    nano-strain/a    ',
     . 'nano-strain/a    nano-strain/a    nano-strain/a    ',
     . 'nano-strain/a        degree       nano-strain/a    km',
     . '        /         /       /')
91    format(f8.3,f9.3,1x,2(f7.1,f6.1),f7.2,f8.1,f7.1,6(f9.1,f8.1),2f9.1,
     .      (f9.1,f8.1),f9.1,2f9.2,i6)
92    format(f8.3,f9.3,2x,7hindex =,i2)

      close(9)
      stop
      end

      subroutine visr_core(yy,xx,min_tau,max_tau,in_tau,
     .                     is_wght,id_wght,wt0,id,indxx)
      implicit real*8 (a-h,l,o-z)
      integer*8 nstn,nslct,ncrp
      parameter (maxsit=5000,maxc=10,maxstn=2000)
      dimension lon(maxsit),lat(maxsit),uxl(maxsit),sxl(maxsit),
     .          uyl(maxsit),syl(maxsit),cxy(maxsit),area(maxsit)
      dimension alon(maxc),alat(maxc),blon(maxc),blat(maxc),
     .          dcs(maxc),dsn(maxc)
      dimension crp_strk1(maxc),crp_strk2(maxc),is0(maxc),icrp(maxc)
      dimension istn(maxstn),azi(maxstn),wght(maxstn)
      dimension a(2*maxstn,6),b(2*maxstn),aa(6,6),aasv(6,6),
     .          ainv(6,6),indx(6),bb(6),v(6)
      common /veldat/nstn,lon,lat,uxl,sxl,uyl,syl,cxy,area
      common /crpdat/ncrp,alon,alat,blon,blat,dcs,dsn
      common /inv_v/bb,ainv,rtau,wtt,nslct,chisq
      common /prmtr/ pi,cov,cutoff_dis,wt_az
*
*  compute relation between interpolation point and creep fault
*
      do i = 1, ncrp

        icrp(i) = 0

        dx = alon(i) - xx
        dy = alat(i) - yy
        ang1 = datan2(dy,dx)/cov
        y1 = -1.0d0*(dy*dcs(i)-dx*dsn(i))

        dx = blon(i) - xx
        dy = blat(i) - yy
        ang2 = datan2(dy,dx)/cov

        if ( ang1.eq.ang2 ) then
          is0(i) = 0                                                         !  interpolation point is located on the delay line of creep fault ?
          goto 20
        endif 
        if ( y1.eq.0.0d0 ) then
          indxx = 4                                                          !  interpolation point is located on crrep fault
          goto 60  
        endif

        is0(i) = y1/abs(y1)

20      if (ang2 .gt. ang1) then
          crp_strk1(i) = ang1
          crp_strk2(i) = ang2
        else
          crp_strk1(i) = ang2
          crp_strk2(i) = ang1
        endif

        dang = crp_strk2(i) - crp_strk1(i)
        if ( dang.gt.180.0d0 ) then
          ang_tmp = crp_strk1(i) + 360.0d0
          crp_strk1(i) = crp_strk2(i)
          crp_strk2(i) = ang_tmp
          icrp(i) = 1
        endif

      enddo
*
*  start data selection.
*
      do 30 itau = min_tau, max_tau, in_tau

        if (itau.eq.max_tau) then
            write(*,*) 'Reached max D but did not find optimal!'
            goto 60
        endif

        indxx = 0
        nslct = 0
        rtau = itau
*       write(*,2000) 'Cutoff dis for D=', rtau, 'is', cutoff_dis*rtau
2000    format(a20,F10.4,a4,F10.4)

        do 10 i = 1, nstn
          if ( id.eq.1.and.lon(i).eq.xx.and.lat(i).eq.yy ) goto 10

          dx = lon(i) - xx
          dy = lat(i) - yy
          dr = dsqrt(dy**2+dx**2)
          ang_s1 = datan2(dy,dx)/cov
          rt = dr/rtau

          if (rt .le. cutoff_dis) then
            do j = 1, ncrp
              if (icrp(j).eq.1.and.ang_s1.lt.0.d0) ang_s1=ang_s1+360.d0
 
              dx = alon(j) - lon(i)
              dy = alat(j) - lat(i)
              y1 = -1.0d0*(dy*dcs(j)-dx*dsn(j))
              
              if ( y1.eq.0.0d0) then
                 is1 = 0
              else
                 is1 = y1/abs(y1)
              endif
             
              if (is1.ne.is0(j)) then
              if (ang_s1.ge.crp_strk1(j).and.
     .          ang_s1.le.crp_strk2(j)) then
                goto 10
              endif
              endif

            enddo
            nslct = nslct + 1
            istn(nslct) = i
          endif
10      continue

      if (nslct.lt.3) then
        indxx = 1
        goto 30
      endif                                                                  ! no interpolation for n_slct < 3

      if (nslct.gt.maxstn) then
        write(6,62) maxstn
        stop
      endif

62    format('# of selected stations > maxstn',i5,
     .   ', please increase maxstn ( visr_core )...')
*
*  estimate azimuthal coverage.
*
      do i = 1, nslct
        dx = lon(istn(i)) - xx
        dy = lat(istn(i)) - yy
        azi(i) = datan2(dy,dx)/cov
      enddo
      azimax = azi(1)
      azimin = azi(1)
      do i = 2, nslct
        azimax = max(azimax,azi(i))
        azimin = min(azimin,azi(i))
      enddo
      dazim1 = azimax - azimin
      if (dazim1.le.180.0d0) then
        indxx = 2
        goto 30
      endif                                                                  ! no interpolation for dazi < pi

      do i = 1, nslct
        if (azi(i).lt.0.d0 ) azi(i) = azi(i) + 360.0d0
        if ( i.eq.1 ) then
          azimax = azi(i)
          azimin = azi(i)
        else
          azimax = max(azimax,azi(i))
          azimin = min(azimin,azi(i))
        endif
      enddo
      dazim2 = azimax - azimin
      if (dazim2.le.180.0d0) then
        indxx = 2
        goto 30
      endif                                                                  ! no interpolation for dazi < pi
*
*  compute data spatial density weighting 
*
      if ( id_wght.eq.1 ) then 
     
        azi_avrg = wt_az*360.0d0/nslct                                       
        azi_tot = (1.0d0+wt_az)*360.0d0                                      

        do i = 1, nslct
          daz1 = 180.0d0
          daz2 = -180.0d0
          do j = 1, nslct
            if (i.eq.j) goto 22
            dazi = azi(j) - azi(i)
            if (dazi.gt.180.0d0) dazi = dazi - 360.0d0
            if (dazi.lt.-180.0d0) dazi = dazi + 360.0d0
            if (dazi.gt.0.0d0.and.dazi.lt.daz1) daz1 = dazi
            if (dazi.lt.0.0d0.and.dazi.gt.daz2) daz2 = dazi
22        continue
          enddo
          wght(i)=(0.5d0*(daz1-daz2)+azi_avrg)*nslct/azi_tot                 ! compute azimuthal weighting
        enddo

      endif

      if ( id_wght.eq.2 ) then

        sum_area = 0.0d0
        do i = 1, nslct
          sum_area = sum_area + area(istn(i))
        enddo
        sum_area = sum_area/nslct

        do i = 1, nslct
          wght(i) = area(istn(i))/sum_area
        enddo

      endif
                  
      npr = 2*nslct
      wtt = 0.0d0
      do i = 1, nslct
        dy = lat(istn(i))-yy
        dx = lon(istn(i))-xx
        dr = dsqrt(dx**2+dy**2)
        if (is_wght .eq. 1) wti = dexp(-(dr/rtau)**2)*wght(i)
        if (is_wght .eq. 2) wti = wght(i)/((dr/rtau)**2+1)
        wtt = wtt + 2*wti

        a(2*i-1,1)=sqrt(wti)/sxl(istn(i))
        a(2*i-1,2)=0.0d0
        a(2*i-1,3)=sqrt(wti)*dx/sxl(istn(i))
        a(2*i-1,4)=sqrt(wti)*dy/sxl(istn(i))
        a(2*i-1,5)=0.0d0
        a(2*i-1,6)=sqrt(wti)*dy/sxl(istn(i))
        a(2*i,1)=0.0d0
        a(2*i,2)=sqrt(wti)/syl(istn(i))
        a(2*i,3)=0.0d0
        a(2*i,4)=sqrt(wti)*dx/syl(istn(i))
        a(2*i,5)=sqrt(wti)*dy/syl(istn(i))
        a(2*i,6)=-1.0*sqrt(wti)*dx/syl(istn(i))

        b(2*i-1)=sqrt(wti)*uxl(istn(i))/sxl(istn(i))
        b(2*i)=sqrt(wti)*uyl(istn(i))/syl(istn(i))
      enddo

      if (wtt .lt. wt0) then
        indxx = 3
        goto 30
      endif                                                                  ! no interpolation for total weight < weighting threshold

      ww = 0.0d0 
      do i = 1, npr
         ww = ww + b(i)**2
      enddo

      do i=1,6
         do j=1,i
            aa(i,j)=0.0d0
            do k=1,npr
               aa(i,j)=aa(i,j)+a(k,i)*a(k,j)
            enddo
            aa(j,i)=aa(i,j)
         enddo
      enddo
         
      do i = 1, 6
         do j = 1, 6
            aasv(i,j) = aa(i,j)
         enddo
      enddo

      do i=1,6
         bb(i)=0.0d0
         do j=1,npr
            bb(i)=bb(i)+a(j,i)*b(j)
         enddo
      enddo
*
*  solve b=n*x.
*
      call dludcmp(aa,6,6,indx,d)
      call dlubksb(aa,6,6,indx,bb)
*
*  invert n.
*
      do i=1,6
         do j=1,6
            ainv(i,j)=0.0d0
         enddo
         ainv(i,i)=1.0d0
      enddo
      do i=1,6
         call dlubksb(aa,6,6,indx,ainv(1,i))
      enddo

      call mat_aba(bb,aasv,6,xaax)
      chisq = ww - xaax
      goto 60

30    continue

60    return
      end

      subroutine cmp_strain(x,covx,emax,emaxsd,emin,eminsd,
     .           taumax,tausd,dexazim,azsd,dilat,ddilat)
      implicit real*8 (a-h,o-z)
      dimension x(3),covx(3,3),v(3)
      common /prmtr/ pi,cov,cutoff_dis,wt_az
*
*  estimate principle strain rates emax, emin, maximum shear tau_max, 
*  and dextral tau_max azimuth
*
      emean=(x(1)+x(3))/2.0d0
      ediff=(x(1)-x(3))/2.0d0
      taumax=dsqrt(x(2)**2+ediff**2)
      emax=emean+taumax
      emin=emean-taumax
      azim=-datan2(x(2),ediff)/cov/2.0d0
      azim=90.0d0+azim
      dexazim=azim+45.0d0-180.0d0
*
*  estimate sigma of tau_max
*
      v(1)=(x(1)-x(3))/4.0d0/taumax
      v(2)=x(2)/taumax
      v(3)=-v(1)
      call mat_aba(v,covx,3,tausd) 
      tausd=dsqrt(tausd)
*
*  estimate sigma of emax
*
      v(1)=0.5d0*(1+(x(1)-x(3))/2.d0/taumax)
      v(2)=x(2)/taumax
      v(3)=0.5d0*(1-(x(1)-x(3))/2.d0/taumax)
      call mat_aba( v,covx,3,emaxsd) 
      emaxsd=dsqrt(emaxsd)
*
*  estimate sigma of emin
*
      v(1)=0.5d0*(1-(x(1)-x(3))/2.0d0/taumax)
      v(2)=-x(2)/taumax
      v(3)=0.5d0*(1+(x(1)-x(3))/2.0d0/taumax)
      call mat_aba(v,covx,3,eminsd) 
      eminsd=dsqrt(eminsd)
*
*  estimate sigma of azimuth
*
      cf=1.0d0/((x(1)-x(3))**2+4.0d0*x(2)**2)
      v(1)=cf*x(2)
      v(2)=-cf*(x(1)-x(3))
      v(3)=-v(1)
      call mat_aba(v,covx,3,azsd) 
      azsd = dsqrt(azsd)/cov
*
*  compute dilatation and its sigma 
*
      dilat=x(1)+x(3)
      ddilat=dsqrt(covx(1,1)+covx(3,3)+2*covx(1,3))

      return
      end

      subroutine mat_aba(a,b,n,aba)
      integer*4 n
      real*8 a(n),b(n,n),aba

      aba=0.0d0
      do i=1,n
         do j=1,n
            aba=aba+a(i)*b(i,j)*a(j)
         enddo
      enddo
      return
      end

      subroutine llxy(slatm,slonm,slat,slon,m)
*
*   subroutine to transfer latitude and longitude into local x y.
*
*   input:
*      slatm: latitude of reference point for local coordinates, in degree.
*      slonm: longitude of reference point for local coordinates, in degree.
*      slat : array coordinate of latitude to be transferred, in degree.
*      slon : array coordinate of longitude to be transferred, in degree.
*      m    : number of points to be transferred.
*   output:
*      slat : y array after transformation in km.
*      slon : x array after transformation in km.
*
      implicit real*8 (a-h,o-z)
      integer*8 m
      dimension t(3,3),vp(3),v(3)
      dimension slat(m),slon(m)
*
      pi = 4.0d0*datan(1.0d0)
      rlatc = slatm*pi/180.0d0
      rlonc = slonm*pi/180.0d0
*
*     calculate local radius of curvature (r) using reference earth nad 
*     1983 (semi-major axis is 6378.137,flattening factor is 1/298.2572)
*
      flat = 1.0d0/298.2572d0
      esq = 2.0d0*flat - flat**2
      q = 1.0d0 - esq*dsin(rlatc)*dsin(rlatc)
      r = 6378.137d0*sqrt(1.0d0 - esq)/q
*
*     construct transformation matrix
*
      t(1,1) =  dsin(rlatc)*dcos(rlonc)
      t(1,2) =  dsin(rlatc)*dsin(rlonc)
      t(1,3) = -dcos(rlatc)
      t(2,1) = -dsin(rlonc)
      t(2,2) =  dcos(rlonc)
      t(2,3) =  0.0d0
      t(3,1) =  dcos(rlatc)*dcos(rlonc)
      t(3,2) =  dcos(rlatc)*dsin(rlonc)
      t(3,3) =  dsin(rlatc)
*
*     calculate xl,yl,dis (vector in local coordinate)
*
      do 9 i = 1,m
        slat(i) = slat(i)*pi/180.0d0
        slon(i) = slon(i)*pi/180.0d0
        v(1) = dcos(slat(i))*dcos(slon(i))
        v(2) = dcos(slat(i))*dsin(slon(i))
        v(3) = dsin(slat(i))
        do j = 1,2
          vp(j) = 0.0d0
          do k = 1,3
            vp(j) = vp(j) + t(j,k)*v(k)
          end do
        end do
        slon(i) =  r*vp(2)
        slat(i) = -r*vp(1)
9     continue

      return
      end

      subroutine dlubksb(a,n,np,indx,b)
      implicit real*8 (a-h,o-z)
      dimension a(np,np),indx(np),b(np)
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.0d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
13        continue
        endif
        b(i)=sum/a(i,i)
14    continue
      return
      end
 
      subroutine dludcmp(a,n,np,indx,d)
      implicit real*8 (a-h,o-z)
      parameter (nmax=1000,tiny=1.0e-20)
      dimension a(np,np),indx(np),vv(nmax)
      d=1.0d0
      do 12 i=1,n
        aamax=0.0d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.0d0) write (6,100)
        vv(i)=1.0d0/aamax
12    continue
100   format('singular matrix.')

      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.0d0
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.0d0)a(j,j)=tiny
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.0.0d0)a(n,n)=tiny
      return
      end

      subroutine cmp_area1(nstn,pos,np,cfa1,ith,area_ith)
      implicit real*8 (a-h,l,o-z)
      integer*8 nstn
      dimension pos(2,nstn),dsmin(np)

      pi = 4.0d0*datan(1.0d0)

      is = 1
      dsmax = 0.0d0
90    dsmin(is) = 1000000.0d0
      do i = 1,nstn
        if ( i.ne.ith) then
          dx = pos(1,ith)-pos(1,i)
          dy = pos(2,ith)-pos(2,i)
          ds = sqrt(dx*dx+dy*dy)
          if ( ds.lt.dsmin(is).and.ds.gt.dsmax ) then
            dsmin(is) = ds
          endif
        endif
      enddo
      dsmax = dsmin(is)

      if ( is.eq.np ) goto 100
      is = is + 1
      goto 90

100   area_ith = 0.0d0

      do i = 1,is
        area_ith = area_ith + dsmin(i)
      enddo
      area_ith = cfa1*area_ith/is/2.0d0
      area_ith = pi*area_ith**2

      return
      end
