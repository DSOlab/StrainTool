subroutine get_voronoi_area_version_init(g_num,g_xy ,voronoi_area)
! Parameters:
!    Input, g_num  integer ( kind = 4 ), the number of  points.
!           g_xy(2,g_num)  real ( kind = 8 ) , the coordinates of input points.
!    output, voronoi_area(g_num) real ( kind = 8 ) , the area of each voronoi cell .
!            the area of the nodes that form the convex hull  is flaged by -1.00

     
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) G_DEGREE(G_NUM), the degree of each Voronoi 
!    cell.
!
!    Local, integer ( kind = 4 ) G_FACE(6*G_NUM), the sequence of vertices to 
!    be used in a traversal of the boundary of the cell associated with each 
!    point.
!
!    Local, integer ( kind = 4 ) G_START(G_NUM), the index in G_FACE of the 
!    first vertex at which to begin a traversal of the boundary of the 
!    cell associated with each point.
!
!    Local, integer ( kind = 4 ) I_NUM, the number of vertices at infinity of the 
!    Voronoi diagram.
!
!    Local, real ( kind = 8 ) I_XY(2,I_NUM), the direction of the
!    vertices at infinity.
!
!    Local, integer ( kind = 4 ) V_NUM, the number of vertices of the Voronoi 
!    diagram.
!    
!    Local , integer ( kind = 4 ) HULL_NUM, the number of nodes that lie on 
!    the convex hull.
!
!    Local , integer ( kind = 4 ) HULL(G_NUM).  Entries 1 through HULL_NUM 
!    contain the indices of the nodes that form the convex hull, in order.
!
 
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  real ( kind = 8 )  g_xy(dim_num,g_num)
  real ( kind = 8 )  voronoi_area(g_num)
  integer ( kind = 4 ) g_num
  
  integer ( kind = 4 ) hull_num
  integer ( kind = 4 ) i_num
  integer ( kind = 4 ) n_num
  integer ( kind = 4 ) v_num
  integer ( kind = 4 ) tri_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k  
  real    ( kind = 8 ) area_temp
  integer ( kind = 4 ) n_temp
  
  integer ( kind = 4 ), allocatable, dimension ( : ) :: g_degree
  integer ( kind = 4 ), allocatable, dimension ( : ) :: g_face
  integer ( kind = 4 ), allocatable, dimension ( : ) :: g_start
  integer ( kind = 4 ), allocatable, dimension ( : ) :: hull
  integer ( kind = 4 ), allocatable, dimension ( :,: ) :: nodtri 
  integer ( kind = 4 ), allocatable, dimension ( :,: ) :: tnbr 
  integer ( kind = 4 ), allocatable, dimension ( : ) ::  indx

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v_xy
  real ( kind = 8 ), allocatable, dimension ( : ) :: cell_x
  real ( kind = 8 ), allocatable, dimension ( : ) :: cell_y
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: i_xy
  
  allocate ( hull(g_num) )
  allocate ( cell_x(20) )
  allocate ( cell_y(20) )
  allocate ( nodtri(3,2*g_num) )
  allocate ( tnbr(3,2*g_num) )
  allocate ( indx(g_num) )
  
  call dtris2 ( g_num, g_xy, tri_num, nodtri, tnbr )  
  call points_hull( g_num,tri_num, nodtri, hull_num, hull )
  
  
  allocate ( i_xy(dim_num,g_num) )
  allocate ( g_degree(g_num) )
  allocate ( g_face(10*g_num) )
  allocate ( g_start(g_num) )
  allocate ( v_xy(dim_num,2*g_num) )
  
 
! get the  data defining the Voronoi diagram.
  call voronoi_data (g_num, g_xy, g_degree, g_start, g_face, v_num, v_xy, &
    i_num, i_xy )


! calculate the area of ecah voronoi cell
!  open  ( unit = 20, file = 'voronoi_xy.txt' )           ! voronoi cell verticals file:

   do i = 1,g_num
      k = g_start(i)

!     write ( 20, * ) ">> ",i                             ! for output

      do j = 1, g_degree(i)
         cell_x(j) = v_xy(1,g_face(k))
         cell_y(j) = v_xy(2,g_face(k))
!        write ( 20, '(2g25.8)' ) v_xy(1,g_face(k)), v_xy(2,g_face(k))          ! for output
         k = k + 1
      end do
!     write ( 20,'(2g25.8)') v_xy(1,g_face(g_start(i))),v_xy(2,g_face(g_start(i)))                ! for output
      n_temp = g_degree(i)
      call areapg ( n_temp , cell_x, cell_y , area_temp )
      voronoi_area(i) = area_temp
   end do
!  close ( 20 ,status = "KEEP" )

  do i = 1, hull_num
     voronoi_area(hull(i)) = -1.0
  end do

  deallocate ( i_xy )
  deallocate ( g_degree )
  deallocate ( g_face )
  deallocate ( g_start )
  deallocate ( v_xy )
  deallocate ( hull )
  deallocate ( cell_x )
  deallocate ( cell_y )

  return
end

function angle_rad_2d ( p1, p2, p3 )
!*****************************************************************************80
!
!! ANGLE_RAD_2D returns the angle in radians swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
!
!        P1
!        /
!       /    
!      /     
!     /  
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), define the rays
!    P1 - P2 and P3 - P2 which define the angle.
!
!    Output, real ( kind = 8 ) ANGLE_RAD_2D, the angle swept out by the rays,
!    in radians.  0 <= ANGLE_RAD_2D < 2 * PI.  If either ray has zero
!    length, then ANGLE_RAD_2D is set to 0.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) angle_rad_2d
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) p(dim_num)
  real    ( kind = 8 ) p1(dim_num)
  real    ( kind = 8 ) p2(dim_num)
  real    ( kind = 8 ) p3(dim_num)

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( all ( p(1:dim_num) == 0.0D+00)  ) then
    return
  end if

  angle_rad_2d = atan2 ( p(2), p(1) )
 if ( angle_rad_2d < 0.0D+00 ) then
    angle_rad_2d = angle_rad_2d + 2.0D+00 * pi
  end if
return
end
subroutine angle_half_2d ( p1, p2, p3, p4 )
!*****************************************************************************80
!
!! ANGLE_HALF_2D finds half an angle in 2D.
!
!
!  Discussion:
!
!    The original angle is defined by the sequence of points P1, P2 and P3.
!
!    The point P4 is calculated so that:
!
!      (P1,P2,P4) = (P1,P2,P3) / 2
!
!        P1
!        /
!       /   P4
!      /  .  
!     / .
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), points defining the angle. 
!
!    Input, real ( kind = 8 ) P4(2), a point defining the half angle.
  implicit none
  integer ( kind = 4 ), parameter :: dim_num = 2
  real    ( kind = 8 ) p1(dim_num)
  real    ( kind = 8 ) p2(dim_num)
  real    ( kind = 8 ) p3(dim_num)
  real    ( kind = 8 ) p4(dim_num)
  integer    ( kind = 4 ) lr
  integer    ( kind = 4 ) lrline
  lr=lrline(p2(1),p2(2),p1(1),p1(2),p3(1),p3(2),0.0D+00)
  if ( lr == 0 ) then 
       p4(1)=p2(1)+(p1(2)-p2(2))/(sqrt ( sum ( ( p1(1:2) - p2(1:2) )**2 ) ))
       p4(2)=p2(2)-(p1(1)-p2(1))/(sqrt ( sum ( ( p1(1:2) - p2(1:2) )**2 ) ))
  else 
       p4(1:2) = 0.5D+00 * ( &
          ( p1(1:2) - p2(1:2) ) / sqrt ( sum ( ( p1(1:2) - p2(1:2) )**2 ) ) &
        + ( p3(1:2) - p2(1:2) ) / sqrt ( sum ( ( p3(1:2) - p2(1:2) )**2 ) ) )
       p4(1:2) = p2(1:2) + p4(1:2) / sqrt ( sum ( p4(1:2)**2 ) )
  end if 
  return
end

subroutine voronoi_data ( g_num, g_xy, g_degree, g_start, g_face, v_num, &
  v_xy, i_num, i_xy )

!*****************************************************************************80
!
!! VORONOI_DATA returns data defining the Voronoi diagram.
!
!  Discussion:
!
!    The routine first determines the Delaunay triangulation.
!
!    The Voronoi diagram is then determined from this information.
!
!    In particular, the circumcenter of each Delaunay triangle
!    is a vertex of a Voronoi polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) G_NUM, the number of generators.
!
!    Input, real ( kind = 8 ) G_XY(2,G_NUM), the point coordinates.
!
!    Output, integer ( kind = 4 ) G_DEGREE(G_NUM), the degree of each 
!    Voronoi cell.
!
!    Output, integer ( kind = 4 ) G_START(G_NUM), the index in G_FACE of the 
!    first vertex at which to begin a traversal of the boundary of the 
!    cell associated with each point.
!
!    Output, integer ( kind = 4 ) G_FACE(6*G_NUM), the sequence of vertices to 
!    be used in a traversal of the boundary of the cell associated with each 
!    point.
!
!    Output, integer ( kind = 4 ) V_NUM, the number of vertices of the Voronoi 
!    diagram.
!
!    Output, real ( kind = 8 ) V_XY(2,V_NUM), the coordinates of the vertices
!    Output, real ( kind = 8 ) I_XY(2,I_NUM), the direction of the
!    vertices at infinity.
!
  implicit none

  integer ( kind = 4 ) g_num
  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) area
  integer ( kind = 4 ) count
  logical, parameter :: debug = .true.
  integer ( kind = 4 ) g
  integer ( kind = 4 ) g_degree(g_num)
  integer ( kind = 4 ) g_face(10*g_num)
  integer ( kind = 4 ) g_next
  integer ( kind = 4 ) g_start(g_num)
  real    ( kind = 8 ) g_xy(dim_num,g_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_num
  integer ( kind = 4 ) i4_wrap
  real    ( kind = 8 ) i_xy(dim_num,g_num)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ix1
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nodtri(3,2*g_num)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) sp1
  integer ( kind = 4 ) s_next
  real    ( kind = 8 ) t(dim_num,3)
  integer ( kind = 4 ) tnbr(3,2*g_num)
  integer ( kind = 4 ) v
  integer ( kind = 4 ) v_inf
  integer ( kind = 4 ) v_next
  integer ( kind = 4 ) v_num
  integer ( kind = 4 ) v_old
  integer ( kind = 4 ) v_save
  real    ( kind = 8 ) v_xy(dim_num,2*g_num)
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
!
!  Compute the Delaunay triangulation.
!
!  do v=1,g_num
!     write ( * ,'(2f14.6)' ) g_xy(1,v) , g_xy(2,v)
!  end do 

  call dtris2 ( g_num, g_xy, v_num, nodtri, tnbr )
!
!  Compute and print the areas of the finite triangles.
!
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) '  Triangle    Area'
!  write ( *, '(a)' ) ' '
!
  do v = 1, v_num
    i1 = nodtri(1,v)
    i2 = nodtri(2,v)
    i3 = nodtri(3,v)

    t(1:dim_num,1:3) = reshape ( (/ &
      g_xy(1:2,i1), g_xy(1:2,i2), g_xy(1:2,i3) /), (/ dim_num, 3 /) )

    call triangle_area_2d ( t, area )

!   write ( *, '(2x,i8,2x,g14.6)' ) v, area

  end do
!
!  Extend the NODTRI data structure, adding fictitious vertices at infinity,
!  so that the Delaunay triangulation can be regarded as covering the
!  entire plane.
!
  call tri_augment ( v_num, nodtri, v_inf )

  if ( debug ) then

!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) '  The generators that form each Delaunay triangle:'
!    write ( *, '(a)' ) '  (Negative values are fictitious nodes at infinity.)'
!    write ( *, '(a)' ) ' '

!    call i4mat_transpose_print ( 3, v_num+v_inf, nodtri, '  Triangle nodes:' )

  end if
!
!  Negative entries in TNBR indicate a semi-infinite Voronoi side.
!  However, DTRIS2 uses a peculiar numbering.  Renumber them.
!
  i_num = 0
  do v = 1, v_num
    do i = 1, 3
    if ( tnbr(i,v) < 0 ) then
        i_num = i_num + 1
        tnbr(i,v) = -i_num
      end if
    end do
  end do

!  if ( debug ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) '  Neighboring triangles of each Delaunay triangle:'
!    write ( *, '(a)' ) '  Negative values indicate no finite neighbor.'
!    write ( *, '(a)' ) ' '

!    call i4mat_transpose_print ( 3, v_num, tnbr, '  Neighbor triangles:' )
!
!  end if
!
!  Determine the degree of each cell.
!
  g_degree(1:g_num) = 0

  do j = 1, v_num + v_inf
    do i = 1, 3
      k = nodtri(i,j)
      if ( 0 < k ) then
        g_degree(k) = g_degree(k) + 1
      end if
    end do
  end do

!  call i4vec_print ( g_num, g_degree, '  Voronoi cell degrees' )
!
!  Each (finite) Delaunay triangle contains a vertex of the Voronoi polygon,
!  at the triangle's circumcenter.
!
  do v = 1, v_num

    i1 = nodtri(1,v)
    i2 = nodtri(2,v)
    i3 = nodtri(3,v)

    t(1:dim_num,1:3) = reshape ( (/ &
      g_xy(1:2,i1), g_xy(1:2,i2), g_xy(1:2,i3) /), (/ dim_num, 3 /) )

    call triangle_circumcenter_2d ( t, v_xy(1:2,v) )

  end do

!  call r8mat_transpose_print ( dim_num, v_num, v_xy, '  The Voronoi vertices:' )
!
!  For each generator G:
!    Determine if its region is infinite.
!      Find a Delaunay triangle containing G.
!      Seek another triangle containing the next node in that triangle.
!
  count = 0
  g_start(1:g_num) = 0

  do g = 1, g_num

    v_next = 0
    do v = 1, v_num + v_inf
      do s = 1, 3
        if ( nodtri(s,v) == g ) then
          v_next = v
          s_next = s
          exit
        end if
      end do
      if ( v_next /= 0 ) then
        exit
      end if
    end do

    v_save = v_next

    do

      s_next = i4_wrap ( s_next + 1, 1, 3 )
      g_next = nodtri(s_next,v_next)

      if ( g_next == g ) then
        s_next = i4_wrap ( s_next + 1, 1, 3 )
        g_next = nodtri(s_next,v_next)
      end if

      v_old = v_next
      v_next = 0
      do v = 1, v_num + v_inf

        if ( v == v_old ) then
          cycle
        end if

        do s = 1, 3

          if ( nodtri(s,v) == g ) then
            sp1 = i4_wrap ( s + 1, 1, 3 )
            if ( nodtri(sp1,v) == g_next ) then
              v_next = v
              s_next = sp1
              exit
            end if
            sp1 = i4_wrap ( s + 2, 1, 3 )
            if ( nodtri(sp1,v) == g_next ) then
              v_next = v
              s_next = sp1
              exit
            end if
          end if

        end do
        if ( v_next /= 0 ) then
          exit
        end if
      end do

      if ( v_next == v_save ) then
        exit
      end if

      if ( v_next == 0 ) then
        v_next = v_old
        exit
      end if

    end do
!
!  Now, starting in the current triangle, V_NEXT, cycle again,
!  and copy the list of nodes into the array.
!
    v_save = v_next

    count = count + 1
    g_start(g) = count
    g_face(count) = v_next

    do
      s_next = i4_wrap ( s_next + 1, 1, 3 )
      g_next = nodtri(s_next,v_next)

      if ( g_next == g ) then
        s_next = i4_wrap ( s_next + 1, 1, 3 )
        g_next = nodtri(s_next,v_next)
      end if

      v_old = v_next
      v_next = 0
      do v = 1, v_num + v_inf
        if ( v == v_old ) then
          cycle
        end if

        do s = 1, 3

          if ( nodtri(s,v) == g ) then
            sp1 = i4_wrap ( s + 1, 1, 3 )
            if ( nodtri(sp1,v) == g_next ) then
              v_next = v
              s_next = sp1
              exit
            end if
            sp1 = i4_wrap ( s + 2, 1, 3 )
            if ( nodtri(sp1,v) == g_next ) then
              v_next = v
              s_next = sp1
              exit
            end if
          end if

        end do
        if ( v_next /= 0 ) then
          exit
        end if
      end do

      if ( v_next == v_save ) then
        exit
      end if

      if ( v_next == 0 ) then
        exit
      end if

      count = count + 1
      g_face(count) = v_next

    end do
  end do
!
!  Mark all the vertices at infinity with a negative sign,
!  so that the data in G_FACE is easier to interpret.
!
  do i = 1, count
    if ( v_num < g_face(i) ) then
      g_face(i) = -g_face(i)
    end if
  end do
!
!  For each (finite) Delaunay triangle, I
!  For each side J,
!
  do i = 1, v_num
    do j = 1, 3
      k = tnbr(j,i)

!  If there is no neighboring triangle on that side,
!  extend a line from the circumcenter of I in the direction of the
!  outward normal to that side.  This is an infinite edge of
!  an infinite Voronoi polygon.
!
      if ( k < 0 ) then

        ix1 = nodtri(j,i)
        x1 = g_xy(1,ix1)
        y1 = g_xy(2,ix1)

        jp1 = i4_wrap ( j+1, 1, 3 )

        ix2 = nodtri(jp1,i)
        x2 = g_xy(1,ix2)
        y2 = g_xy(2,ix2)
!
!  Compute the direction I_XY(1:2,-K).
!
!        call line_exp_normal_2d ( g_xy(1:2,ix1), g_xy(1:2,ix2), i_xy(1:2,-k) )

      end if

    end do
  end do

  return
end

subroutine dtris2 ( point_num, point_xy, tri_num, tri_vert, tri_nabe )
!*****************************************************************************80
!
!! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!  Discussion:
!
!    The routine constructs the Delaunay triangulation of a set of 2D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2004
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of vertices.
!
!    Input/output, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the vertices.
!    On output, the vertices have been sorted into dictionary order.
!
!    Output, integer ( kind = 4 ) TRI_NUM, the number of triangles in the 
!    triangulation;  TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the 
!    number of boundary vertices.
!
!    Output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the nodes that make up 
!    each triangle.  The elements are indices of POINT_XY.  The vertices of the 
!    triangles are in counter clockwise order.
!
!
!    Output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbor 
!    list.  Positive elements are indices of TIL; negative elements are used 
!    for links of a counter clockwise linked list of boundary edges; 
!    LINK = -(3*I + J-1) where I, J = triangle, edge index; TRI_NABE(J,I) refers
!    to the neighbor along edge from vertex J to J+1 (mod 3).
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indx(point_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  real ( kind = 8 ) point_xy(dim_num,point_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) stack(point_num)
  integer ( kind = 4 ) t
  real ( kind = 8 ) tol
  integer ( kind = 4 ) top
  integer ( kind = 4 ) tri_nabe(3,point_num*2)
  integer ( kind = 4 ) tri_num
  integer ( kind = 4 ) tri_vert(3,point_num*2)

  tol = 100.0D+00 * epsilon ( tol )

  ierr = 0
!
!  Sort the vertices by increasing (x,y).
!
  call r82vec_sort_heap_index_a ( point_num, point_xy, indx )

  call r82vec_permute ( point_num, point_xy, indx )
!
!  Make sure that the data points are "reasonably" distinct.
!
  m1 = 1

  do i = 2, point_num

    m = m1
    m1 = i

    k = 0

    do j = 1, 2

      cmax = max ( abs ( point_xy(j,m) ), abs ( point_xy(j,m1) ) )

      if ( tol * ( cmax + 1.0D+00 ) &
           < abs ( point_xy(j,m) - point_xy(j,m1) ) ) then
        k = j
        exit
      end if

    end do

    if ( k == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a,i8)' ) '  Fails for point number I = ', i
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  M1 = ', m1
      write ( *, '(a,2g14.6)' ) '  X,Y(M)  = ', point_xy(1,m), point_xy(2,m)
      write ( *, '(a,2g14.6)' ) '  X,Y(M1) = ', point_xy(1,m1), point_xy(2,m1)
      ierr = 224
      return
    end if

  end do
!
!  Starting from points M1 and M2, search for a third point M that
!  makes a "healthy" triangle (M1,M2,M)
!
  m1 = 1
  m2 = 2
  j = 3

  do

    if ( point_num < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      ierr = 225
      return
    end if

    m = j

    lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
      point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1

  end do
!
!  Set up the triangle information for (M1,M2,M), and for any other
!  triangles you created because points were collinear with M1, M2.
!
  tri_num = j - 2

  if ( lr == -1 ) then

    tri_vert(1,1) = m1
    tri_vert(2,1) = m2
    tri_vert(3,1) = m
    tri_nabe(3,1) = -3

    do i = 2, tri_num

      m1 = m2
      m2 = i+1
      tri_vert(1,i) = m1
      tri_vert(2,i) = m2
      tri_vert(3,i) = m
      tri_nabe(1,i-1) = -3 * i
      tri_nabe(2,i-1) = i
      tri_nabe(3,i) = i - 1
    end do

    tri_nabe(1,tri_num) = -3 * tri_num - 1
    tri_nabe(2,tri_num) = -5
    ledg = 2
    ltri = tri_num

  else

    tri_vert(1,1) = m2
    tri_vert(2,1) = m1
    tri_vert(3,1) = m
    tri_nabe(1,1) = -4

    do i = 2, tri_num
      m1 = m2
      m2 = i+1
      tri_vert(1,i) = m2
      tri_vert(2,i) = m1
      tri_vert(3,i) = m
      tri_nabe(3,i-1) = i
      tri_nabe(1,i) = -3 * i - 3
      tri_nabe(2,i) = i - 1
    end do

    tri_nabe(3,tri_num) = -3 * tri_num
    tri_nabe(2,1) = -3 * tri_num - 2
    ledg = 2
    ltri = 1

  end if
!
!
!  Insert the vertices one at a time from outside the convex hull,
!  determine visible boundary edges, and apply diagonal edge swaps until
!  Delaunay triangulation of vertices (so far) is obtained.
!
  top = 0

  do i = j+1, point_num

    m = i
    m1 = tri_vert(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = tri_vert(ledg+1,ltri)
    else
      m2 = tri_vert(1,ltri)
    end if

    lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
      point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

    if ( 0 < lr ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -tri_nabe(ledg,ltri)
      rtri = l / 3
      redg = mod(l,3) + 1
    end if

    call vbedg ( point_xy(1,m), point_xy(2,m), point_num, point_xy, tri_num, &
      tri_vert, tri_nabe, ltri, ledg, rtri, redg )

    n = tri_num + 1
    l = -tri_nabe(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -tri_nabe(e,t)
      m2 = tri_vert(e,t)

      if ( e <= 2 ) then
        m1 = tri_vert(e+1,t)
      else
        m1 = tri_vert(1,t)
      end if

      tri_num = tri_num + 1
      tri_nabe(e,t) = tri_num
      tri_vert(1,tri_num) = m1
      tri_vert(2,tri_num) = m2
      tri_vert(3,tri_num) = m
      tri_nabe(1,tri_num) = t
      tri_nabe(2,tri_num) = tri_num - 1
      tri_nabe(3,tri_num) = tri_num + 1
      top = top + 1

      if ( point_num < top ) then
        ierr = 8
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
        write ( *, '(a)' ) '  Stack overflow.'
        return
      end if

      stack(top) = tri_num

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    tri_nabe(ledg,ltri) = -3 * n - 1
    tri_nabe(2,n) = -3 * tri_num - 2
    tri_nabe(3,tri_num) = -l
    ltri = n
    ledg = 2

    call swapec ( m, top, ltri, ledg, point_num, point_xy, tri_num, &
      tri_vert, tri_nabe, stack, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a)' ) '  Error return from SWAPEC.'
      return
    end if

  end do
!
!  Now account for the sorting that we did.
!
  do i = 1, 3
    do j = 1, tri_num
      tri_vert(i,j) = indx ( tri_vert(i,j) )
    end do
  end do

  call perm_inv ( point_num, indx )

  call r82vec_permute ( point_num, point_xy, indx )

  return
end

subroutine  areapg ( nvrt, xc, yc, area )
  implicit none
  integer  (  kind = 4 ) nvrt

  real     (  kind = 8 ) area
  integer  (  kind = 4 ) i
  real     (  kind = 8 ) sum2
  real ( kind = 8 ) xc(10)
  real ( kind = 8 ) yc(10)
  sum2 = xc(1) * ( yc(2) - yc(nvrt) )

  do i = 2, nvrt-1
    sum2 = sum2 + xc(i) * ( yc(i+1) - yc(i-1) )
  end do

  sum2 = sum2 + xc(nvrt) * ( yc(1) - yc(nvrt-1) )

  area = ABS(sum2*0.5)
!   write( *,'(a)' ) 'area cal sus~'
  return
end



function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) wide

  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i4_wrap = ilo
  else
    i4_wrap = ilo + i4_modp ( ival - ilo, wide )
  end if

  return
end


function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

!*****************************************************************************80
!
!! LRLINE determines if a point is left of, right or, or on a directed line.
!
!  Discussion:
!
!    The directed line is parallel to, and at a signed distance DV from
!    a directed base line from (XV1,YV1) to (XV2,YV2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XU, YU, the coordinates of the point whose
!    position relative to the directed line is to be determined.
!
!    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, the coordinates of two points
!    that determine the directed base line.
!
!    Input, real ( kind = 8 ) DV, the signed distance of the directed line
!    from the directed base line through the points (XV1,YV1) and (XV2,YV2).
!    DV is positive for a line to the left of the base line.
!
!    Output, integer ( kind = 4 ) LRLINE, the result:
!    +1, the point is to the right of the directed line;
!     0, the point is on the directed line;
!    -1, the point is to the left of the directed line.
!
  implicit none

  real    ( kind = 8 ) dv
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dxu
  real    ( kind = 8 ) dy
  real    ( kind = 8 ) dyu
  integer ( kind = 4 ) lrline
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) tolabs
  real    ( kind = 8 ) xu
  real    ( kind = 8 ) xv1
  real    ( kind = 8 ) xv2
  real    ( kind = 8 ) yu
  real    ( kind = 8 ) yv1
  real    ( kind = 8 ) yv2

  tol = 100.0D+00 * epsilon ( tol )

  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
    abs ( dyu ), abs ( dv ) )

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

  if ( tolabs < t ) then
    lrline = 1
  else if ( -tolabs <= t ) then
    lrline = 0
  else
    lrline = -1
  end if

  return
end

subroutine tri_augment ( v_num, nodtri, v_inf )

!*****************************************************************************80
!
!! TRI_AUGMENT augments the triangle data using vertices at infinity.
!
!  Discussion:
!
!    The algorithm simply looks at the list of triangle edges stored
!    in NODTRI, and determines which edges, of the form (P1,P2), do
!    not have a matching (P2,P1) occurrence.  These correspond to
!    boundary edges of the convex hull.  To simplify our computations,
!    we adjust the NODTRI array to accommodate an extra triangle with 
!    one vertex at infinity for each such unmatched edge.
!
!    The algorithm used here is ruinously inefficient for large V_NUM.
!    Assuming that this data structure modification is the way to go,
!    the routine should be rewritten to determine the boundary edges
!    more efficiently.
!
!    The fictitious vertices at infinity show up in the augmenting
!    rows of the NODTRI array with negative indices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) V_NUM, the number of Voronoi vertices.
!
!    Input/output, integer ( kind = 4 ) NODTRI(3,*), the list of nodes that
!    comprise each Delaunay triangle.  On input, there are V_NUM
!    sets of this data.  On output, for every pair of nodes (P1,P2) 
!    for which the pair (P2,P1) does not occur, an augmenting triangle 
!    has been created with exactly this edge (plus a vertex at infinity).
!    On output, there are V_NUM + V_INF sets of data.
!
!    Output, integer ( kind = 4 ) V_INF, the number of augmenting triangles and 
!    vertices at infinity that were created.
!
  implicit none

  logical              found
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) nodtri(3,*)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) t
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) v
  integer ( kind = 4 ) v_inf
  integer ( kind = 4 ) v_num
  integer ( kind = 4 ) v2

  v_inf = 0

  do v = 1, v_num
    do i = 1, 3

      s = nodtri(i,v)
      ip1 = i4_wrap ( i + 1, 1, 3 )
      t = nodtri(ip1,v)

      found = .false.

      do v2 = 1, v_num

        do i2 = 1, 3
          s2 = nodtri(i2,v2)
          ip1 = i4_wrap ( i2 + 1, 1, 3 )
          t2 = nodtri(ip1,v2)
          if ( s == t2 .and. t == s2 ) then
            found = .true.
            exit
          end if
        end do

        if ( found ) then
          exit
        end if

      end do

      if ( .not. found ) then
        v_inf = v_inf + 1
        nodtri(1:3,v_num+v_inf) = (/ -v_inf, t, s /)
      end if

    end do
  end do

!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) 'TRI_AUGMENT:'
!  write ( *, '(a,i8)' ) '  Number of boundary triangles = ', v_inf

  return
end
subroutine triangle_circumcenter_2d ( t, center )

!*****************************************************************************80
!
!! TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
!
!  Discussion:
!
!    The circumcenter of a triangle is the center of the circumcircle, the
!    circle that passes through the three vertices of the triangle.
!
!    The circumcircle contains the triangle, but it is not necessarily the
!    smallest triangle to do so.
!
!    If all angles of the triangle are no greater than 90 degrees, then
!    the center of the circumscribed circle will lie inside the triangle.
!    Otherwise, the center will lie outside the triangle.
!
!    The circumcenter is the intersection of the perpendicular bisectors
!    of the sides of the triangle.
!
!    In geometry, the circumcenter of a triangle is often symbolized by "O".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) CENTER(2), the circumcenter of the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) asq
  real ( kind = 8 ) bot
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) csq
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ) top(dim_num)

  asq = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2
  csq = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2

  top(1) =    ( t(2,2) - t(2,1) ) * csq - ( t(2,3) - t(2,1) ) * asq
  top(2) =  - ( t(1,2) - t(1,1) ) * csq + ( t(1,3) - t(1,1) ) * asq

  bot  =  ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) ) &
        - ( t(2,3) - t(2,1) ) * ( t(1,2) - t(1,1) )

  center(1:2) = t(1:2,1) + 0.5D+00 * top(1:2) / bot

  return
end

subroutine vbedg ( x, y, point_num, point_xy, tri_num, tri_vert, tri_nabe, &
  ltri, ledg, rtri, redg )
 implicit none

  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  logical              ldone
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  real    ( kind = 8 ) point_xy(2,point_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) t
  integer ( kind = 4 ) tri_nabe(3,tri_num)
  integer ( kind = 4 ) tri_vert(3,tri_num)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
!
!  Find the rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor information.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

  do

    l = -tri_nabe(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = tri_vert(e,t)

    if ( e <= 2 ) then
      b = tri_vert(e+1,t)
    else
      b = tri_vert(1,t)
    end if

    lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
      point_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

    rtri = t
    redg = e

  end do

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = tri_vert(e,t)
    e = i4_wrap ( e-1, 1, 3 )

    do while ( 0 < tri_nabe(e,t) )

      t = tri_nabe(e,t)

      if ( tri_vert(1,t) == b ) then
        e = 3
      else if ( tri_vert(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do
    a = tri_vert(e,t)

    lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
       point_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end
function i4_modp ( i, j )
 implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end

subroutine triangle_area_2d ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_2D computes the area of a triangle in 2D.
!
!  Discussion:
!
!    If the triangle's vertices are given in counter clockwise order,
!    the area will be positive.  If the triangle's vertices are given
!    in clockwise order, the area will be negative!
!
!    An earlier version of this routine always returned the absolute
!    value of the computed area.  I am convinced now that that is
!    a less useful result!  For instance, by returning the signed 
!    area of a triangle, it is possible to easily compute the area 
!    of a nonconvex polygon as the sum of the (possibly negative) 
!    areas of triangles formed by node 1 and successive pairs of vertices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) AREA, the area of the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ) t(dim_num,3)

  area = 0.5D+00 * ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

  return
end

subroutine points_hull( g_num,v_num, nodtri, v_inf, hull )
  implicit none

  logical              found
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) nodtri(3,*)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) t
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) v
  integer ( kind = 4 ) v_inf
  integer ( kind = 4 ) v_num
  integer ( kind = 4 ) v2
  integer ( kind = 4 ) g_num
  integer ( kind = 4 ) hull(g_num)
  logical              point_flag(g_num)
  integer ( kind = 4 ) k
  v_inf = 0
  do i = 1,g_num
     point_flag(i)=.false.
  end do

  do v = 1, v_num
    do i = 1, 3

      s = nodtri(i,v)
      ip1 = i4_wrap ( i + 1, 1, 3 )
      t = nodtri(ip1,v)

      found = .false.

      do v2 = 1, v_num

        do i2 = 1, 3
          s2 = nodtri(i2,v2)
          ip1 = i4_wrap ( i2 + 1, 1, 3 )
          t2 = nodtri(ip1,v2)
          if ( s == t2 .and. t == s2 ) then
            found = .true.
            exit
          end if
        end do

        if ( found ) then
          exit
        end if

      end do

      if ( .not. found ) then
        v_inf = v_inf + 1
!        nodtri(1:3,v_num+v_inf) = (/ -v_inf, t, s /)
         point_flag(t) = .true.
         point_flag(s) = .true.
      end if

    end do
  end do
  k = 0 
  do i = 1,g_num
     if ( point_flag(i) ) then 
       k= k+1
       hull(k) = i
     end if
  end do
    
  return
end


subroutine r82vec_permute ( n, a, p )

!*****************************************************************************80
!
!! R82VEC_PERMUTE permutes an R82VEC in place.
!
!  Discussion:
!
!    This routine permutes an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) a_temp(2)
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp(1:2) = a(1:2,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
          stop
        end if

        if ( iget == istart ) then
          a(1:2,iput) = a_temp(1:2)
          exit
        end if

        a(1:2,iput) = a(1:2,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine r82vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(1:2,INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call R82VEC_PERMUTE ( N, A, INDX )
!
!    after which A(1:2,I), I = 1 to N is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(2,N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(1:2,INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) aval(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then
    indx(1) = 1
    return
  end if

  call i4vec_indicator ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval(1:2) = a(1:2,indxt)

    else

      indxt = indx(ir)
      aval(1:2) = a(1:2,indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
             ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
               a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
          j = j + 1
        end if
      end if

      if (   aval(1) <  a(1,indx(j)) .or. &
           ( aval(1) == a(1,indx(j)) .and. &
             aval(2) <  a(2,indx(j)) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine swapec ( i, top, btri, bedg, point_num, point_xy, tri_num, &
  tri_vert, tri_nabe, stack, ierr )

!*****************************************************************************80
!
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!  Discussion:
!
!    The routine swaps diagonal edges in a 2D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay,
!    given that I is the index of the new vertex added to the triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the new vertex.
!
!    Input/output, integer ( kind = 4 ) TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input/output, integer ( kind = 4 ) BTRI, BEDG; on input, if positive, are 
!    the triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates
!    of the points.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input/output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the triangle 
!    incidence list.  May be updated on output because of swaps.
!
!    Input/output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle 
!    neighbor list; negative values are used for links of the counter-clockwise 
!    linked list of boundary edges;  May be updated on output because of swaps.
!
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer ( kind = 4 ) STACK(MAXST); on input, entries 1 through 
!    TOP contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer IERR is set to 8 for abnormal return.
!
  implicit none

  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bedg
  integer ( kind = 4 ) btri
  integer ( kind = 4 ) c
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) e
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fm1
  integer ( kind = 4 ) fp1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  real    ( kind = 8 ) point_xy(2,point_num)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) stack(point_num)
  integer ( kind = 4 ) swap
  integer ( kind = 4 ) t
  integer ( kind = 4 ) top
  integer ( kind = 4 ) tri_nabe(3,tri_num)
  integer ( kind = 4 ) tri_vert(3,tri_num)
  integer ( kind = 4 ) tt
  integer ( kind = 4 ) u
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  x = point_xy(1,i)
  y = point_xy(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( tri_vert(1,t) == i ) then
      e = 2
      b = tri_vert(3,t)
    else if ( tri_vert(2,t) == i ) then
      e = 3
      b = tri_vert(1,t)
    else
      e = 1
      b = tri_vert(2,t)
    end if

    a = tri_vert(e,t)
    u = tri_nabe(e,t)

    if ( tri_nabe(1,u) == t ) then
      f = 1
      c = tri_vert(3,u)
    else if ( tri_nabe(2,u) == t ) then
      f = 2
      c = tri_vert(1,u)
    else
      f = 3
      c = tri_vert(2,u)
    end if

    swap = diaedg ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,c), &
      point_xy(2,c), point_xy(1,b), point_xy(2,b) )

    if ( swap == 1 ) then

      em1 = i4_wrap ( e - 1, 1, 3 )
      ep1 = i4_wrap ( e + 1, 1, 3 )
      fm1 = i4_wrap ( f - 1, 1, 3 )
      fp1 = i4_wrap ( f + 1, 1, 3 )

      tri_vert(ep1,t) = c
      tri_vert(fp1,u) = i
      r = tri_nabe(ep1,t)
      s = tri_nabe(fp1,u)
      tri_nabe(ep1,t) = u
      tri_nabe(fp1,u) = t
      tri_nabe(e,t) = s
      tri_nabe(f,u) = r

      if ( 0 < tri_nabe(fm1,u) ) then
        top = top + 1
        stack(top) = u
      end if

      if ( 0 < s ) then

        if ( tri_nabe(1,s) == u ) then
          tri_nabe(1,s) = t
        else if ( tri_nabe(2,s) == u ) then
          tri_nabe(2,s) = t
        else
          tri_nabe(3,s) = t
        end if

        top = top + 1

        if ( point_num < top ) then
          ierr = 8
          return
        end if

        stack(top) = t

      else

        if ( u == btri .and. fp1 == bedg ) then
          btri = t
          bedg = e
        end if

        l = - ( 3 * t + e - 1 )
        tt = t
        ee = em1

        do while ( 0 < tri_nabe(ee,tt) )

          tt = tri_nabe(ee,tt)

          if ( tri_vert(1,tt) == a ) then
            ee = 3
          else if ( tri_vert(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tri_nabe(ee,tt) = l

      end if

      if ( 0 < r ) then

        if ( tri_nabe(1,r) == t ) then
          tri_nabe(1,r) = u
        else if ( tri_nabe(2,r) == t ) then
          tri_nabe(2,r) = u
        else
          tri_nabe(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( 0 < tri_nabe(ee,tt) )

          tt = tri_nabe(ee,tt)

          if ( tri_vert(1,tt) == b ) then
            ee = 3
          else if ( tri_vert(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tri_nabe(ee,tt) = l

      end if

    end if

  end do

  return
end
  
function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! DIAEDG chooses a diagonal edge.
!
!  Discussion:
!
!    The routine determines whether 0--2 or 1--3 is the diagonal edge
!    that should be chosen, based on the circumcircle criterion, where
!    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
!    quadrilateral in counterclockwise order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    coordinates of the vertices of a quadrilateral, given in
!    counter clockwise order.
!
!    Output, integer ( kind = 4 ) DIAEDG, chooses a diagonal:
!    +1, if diagonal edge 02 is chosen;
!    -1, if diagonal edge 13 is chosen;
!     0, if the four vertices are cocircular.
!
  implicit none

  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  integer ( kind = 4 ) diaedg
  real ( kind = 8 ) dx10
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx30
  real ( kind = 8 ) dx32
  real ( kind = 8 ) dy10
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy30
  real ( kind = 8 ) dy32
  real ( kind = 8 ) s
  real ( kind = 8 ) tol
  real ( kind = 8 ) tola
  real ( kind = 8 ) tolb
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  tol = 100.0D+00 * epsilon ( tol )

  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2

  tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
  tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

  ca = dx10 * dx30 + dy10 * dy30
  cb = dx12 * dx32 + dy12 * dy32

  if ( tola < ca .and. tolb < cb ) then

    diaedg = -1

  else if ( ca < -tola .and. cb < -tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( tola < s ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end


subroutine perm_inv ( n, p )

!*****************************************************************************80
!
!! PERM_INV inverts a permutation "in place".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard 
!    index form.  On output, P describes the inverse permutation
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) is
  integer ( kind = 4 ) p(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INV - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop
  end if

  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    is = -sign ( 1, p(i) )
    p(i) = sign ( p(i), is )

  end do

  do i = 1, n

    i1 = -p(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = p(i1)
        p(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do

  return
end

subroutine ch_cap ( c )
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
