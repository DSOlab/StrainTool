      subroutine visr_core(
      yy, <- current latitude of cell centre
      xx, <- current longtitude of cell centre
      min_tau,
      max_tau,
      in_tau, <- minimum, maximum, and incremental spatial smoothing constants (km)
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
*  start data selection.
*
      do 30 itau = min_tau, max_tau, in_tau

        if (itau.eq.max_tau) goto 60

        indxx = 0
        nslct = 0
        rtau = itau

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
