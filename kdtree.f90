
Module kd_tree

  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel.  Only the Euclidean metric works so far.


  ! Institute for Nonlinear science, UC San Diego, 1997
  ! kennel@lyapunov.ucsd.edu
  ! .. Parameters ..
  ! you choose this.
  Integer, Parameter :: bucket_size = 10
  ! ..
  ! .. Derived Type Declarations ..
  ! Global information about the tree
  ! pointer to the actual
  ! data array
  ! dimensionality and total # of points
  ! permuted index into the
  ! data, so that
  ! indexes[l..u]
  ! of some bucket represent the indexes of the actual
  ! points in that bucket.
  ! root pointer of the tree
  ! an internal tree node
  ! the dimension to cut
  ! where to cut the dimension
  ! indices of points included in this node,
  ! referring back to indices
  ! child pointers
  ! One of these is created for each search.
  ! best squared distance found so far
  ! best index found so far
  ! Query vector
  ! indexes of best found so far
  ! squared distances found so far
  ! how many best distances we are searching
  ! for
  ! how many have been found so far, i.e. il(1:nfound)
  ! are the only valid indexes.
  ! exclude points within
  ! 'correltime'
  ! of 'centeridx'

  Type :: tree_node
      Integer :: dnum
      Real :: val
      Integer :: l, u
      Type (tree_node), Pointer :: left, right
  End Type tree_node

  Type :: tree_master_record
      Real, Dimension (:,:), Pointer :: the_data
      Integer :: dim, n
      Integer, Dimension (:), Pointer :: indexes
      Type (tree_node), Pointer :: root
  End Type tree_master_record

  Type :: tree_search_record
      Private
      Real :: bsd
      Real, Dimension (:), Pointer :: qv
      Integer, Dimension (:), Pointer :: il
      Real, Dimension (:), Pointer :: dsl
      Integer :: nbst, nfound
      Integer :: centeridx, correltime
      Logical :: linfinity_metric
  End Type tree_search_record
  ! ..
  ! set to true if we're doing linfinity metric
Contains

  Subroutine destroy_tree(tp)
    ! Deallocates all memory for the tree, except input data matrix
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    Call destroy_node(tp%root)
    Deallocate (tp%indexes)
    tp%indexes => NULL() !#Nullify (tp%indexes)
    Return
  end Subroutine destroy_tree
  Recursive Subroutine destroy_node(np)
    ! .. Structure Arguments ..
    Type (tree_node), Pointer :: np
    ! ..
    ! ..
    If (ASSOCIATED(np%left)) Then
       Call destroy_node(np%left)
       Deallocate (np%left)
       np%left => NULL() 
    End If
    If (ASSOCIATED(np%right)) Then
       Call destroy_node(np%right)
       Deallocate (np%right)
       np%right => NULL() 
    End If
    Return
  End Subroutine destroy_node


  Function create_tree(input_data) Result (master_record)
    ! create the actual tree structure, given an input array of data.
    ! Arguments
    ! .. Function Return Value ..
    Type (tree_master_record), Pointer :: master_record
    ! ..
    ! .. Array Arguments ..
    Real, Target :: input_data(:,:)
    ! ..
    ! ..
    Allocate (master_record)
    master_record%the_data => input_data
    ! pointer assignment
    master_record%n = SIZE(input_data,1)
    master_record%dim = SIZE(input_data,2)
    Print *, 'Creating KD tree with N = ', master_record%n, &
     ' and dim = ', master_record%dim
    Call build_tree(master_record)

  Contains

    Subroutine build_tree(tp)
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Local Scalars ..
      Integer :: j
      ! ..
      Allocate (tp%indexes(tp%n))
      Do j = 1, tp%n
         tp%indexes(j) = j
      End Do
      tp%root => build_tree_for_range(tp,1,tp%n)
    End Subroutine build_tree

    Recursive Function build_tree_for_range(tp,l,u) Result (res)
      ! .. Function Return Value ..
      Type (tree_node), Pointer :: res
      ! ..
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer, Intent (In) :: l, u
      ! ..
      ! .. Local Scalars ..
      Integer :: c, m
      ! ..
      If ((u-l)<=bucket_size) Then
         Allocate (res)
         res%dnum = 0
         res%val = 0.0
         res%l = l
         res%u = u
         Nullify (res%left,res%right)
      Else
         Allocate (res)
         c = most_spread_coordinate(tp,l,u)
         m = (l+u)/2
         Call select_on_coordinate(tp,c,m,l,u)
         ! moves indexes around
         res%dnum = c
         res%val = tp%the_data(tp%indexes(m),c)
         res%l = l
         res%u = u
         res%left => build_tree_for_range(tp,l,m)
         res%right => build_tree_for_range(tp,m+1,u)
      End If
    End Function build_tree_for_range

    Subroutine select_on_coordinate(tp,c,k,li,ui)
      ! Move elts of ind around between l and u, so that the kth
      ! element
      ! is >= those below, <= those above, in the coordinate c.
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer, Intent (In) :: c, k, li, ui
      ! ..
      ! .. Local Scalars ..
      Integer :: i, l, m, s, t, u
      ! ..
      ! .. Local Arrays ..
      Real, Pointer :: v(:,:)
      Integer, Pointer :: ind(:)
      ! ..
      v => tp%the_data
      ind => tp%indexes
      l = li
      u = ui
      Do While (l<u)
         t = ind(l)
         m = l
         Do i = l + 1, u
            If (v(ind(i),c)<v(t,c)) Then
               m = m + 1
               s = ind(m)
               ind(m) = ind(i)
               ind(i) = s
            End If
         End Do
         s = ind(l)
         ind(l) = ind(m)
         ind(m) = s
         If (m<=k) l = m + 1
         If (m>=k) u = m - 1
      End Do
    End Subroutine select_on_coordinate

    Function most_spread_coordinate(tp,l,u) Result (res)
      ! Of indices in l..u find the axis which has the largest spread,
      ! and
      ! return its index
      ! .. Function Return Value ..
      Integer :: res
      ! ..
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer, Intent (In) :: l, u
      ! ..
      ! .. Local Scalars ..
      Real :: bsp, sp
      Integer :: i = 0
      ! ..
      res = 0
      bsp = -1.0
      Do i = 1, tp%dim
         sp = spread_in_coordinate(tp,i,l,u)
         If (sp>bsp) Then
            res = i
            bsp = sp
         End If
      End Do
    End Function most_spread_coordinate

    Function spread_in_coordinate(tp,c,l,u) Result (res)
      ! the spread in coordinate 'c', between l and u
      ! for easier local access
      ! ibid
      ! .. Function Return Value ..
      Real :: res
      ! ..
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer, Intent (In) :: c, l, u
      ! ..
      ! .. Local Scalars ..
      Real :: last, lmax, lmin, max, min, t
      Integer :: i
      ! ..
      ! .. Local Arrays ..
      Real, Pointer :: v(:,:)
      Integer, Pointer :: ind(:)
      ! ..
      v => tp%the_data
      ind => tp%indexes
      min = v(ind(l),c)
      max = min
      Do i = l + 2, u, 2
         lmin = v(ind(i-1),c)
         lmax = v(ind(i),c)
         If (lmin>lmax) Then
            t = lmin
            lmin = lmax
            lmax = t
         End If
         If (min>lmin) min = lmin
         If (max<lmax) max = lmax
      End Do
      If (i==u+1) Then
         last = v(ind(u),c)
         If (min>last) min = last
         If (max<last) max = last
      End If
      res = max - min
      ! write *, "Returning spread in coordinate with res = ", res
    End Function spread_in_coordinate
  End Function create_tree

  ! Search routines:

  ! * n_nearest_to(tp, qv, n, indexes, distances)

  ! Find the 'n' vectors in the tree nearest to 'qv' in euclidean norm
  ! returning their indexes and distances in 'indexes' and 'distances'
  ! arrays already allocated passed to the subroutine.
  Subroutine n_nearest_to(tp,qv,n,indexes,distances)
    ! find the 'n' nearest neighbors to 'qv', returning
    ! their indexes in indexes and squared Euclidean distances in
    ! 'distances'
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    ! .. Scalar Arguments ..
    Integer, Intent (In) :: n
    ! ..
    ! .. Array Arguments ..
    Real, Target :: distances(n)
    Real, Target, Intent (In) :: qv(:)
    Integer, Target :: indexes(n)
    ! ..
    ! .. Local Structures ..
    Type (tree_search_record), Pointer :: psr
    Type (tree_search_record), Target :: sr
    ! ..
    ! ..
    ! the largest real number
    sr%bsd = HUGE(1.0)
    sr%qv => qv
    sr%nbst = n
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%linfinity_metric = .False. ! change if you want.
    sr%dsl => distances
    sr%il => indexes
    sr%dsl = HUGE(sr%dsl)    ! set to huge positive values
    sr%il = -1               ! set to invalid indexes
    psr => sr                ! in C this would be psr = &sr
    Call n_search(tp,psr,tp%root)
    Return
  End Subroutine n_nearest_to

  ! Another main routine you call from your user program


  Subroutine n_nearest_to_around_point(tp,idxin,correltime,n,indexes, &
   distances)
    ! find the 'n' nearest neighbors to point 'idxin', returning
    ! their indexes in indexes and squared Euclidean distances in
    ! 'distances'
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    ! .. Scalar Arguments ..
    Integer, Intent (In) :: correltime, idxin, n
    ! ..
    ! .. Array Arguments ..
    Real, Target :: distances(n)
    Integer, Target :: indexes(n)
    ! ..
    ! .. Local Structures ..
    Type (tree_search_record), Pointer :: psr
    Type (tree_search_record), Target :: sr
    ! ..
    ! .. Local Arrays ..
    Real, Allocatable, Target :: qv(:)
    ! ..
    ! ..
    Allocate (qv(tp%dim))
    qv = tp%the_data(idxin,:) ! copy the vector
    sr%bsd = HUGE(1.0)       ! the largest real number
    sr%qv => qv
    sr%nbst = n
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime
    sr%linfinity_metric = .False.
    sr%dsl => distances
    sr%il => indexes
    sr%dsl = HUGE(sr%dsl)    ! set to huge positive values
    sr%il = -1               ! set to invalid indexes
    psr => sr                ! in C this would be psr = &sr
    Call n_search(tp,psr,tp%root)
    Deallocate (qv)
    Return
  End Subroutine n_nearest_to_around_point

  Function bounded_square_distance(tp,qv,ii,bound) Result (res)
    ! distance between v[i,*] and qv[*] if less than or equal to bound,
    ! -1.0 if greater than than bound.
    ! .. Function Return Value ..
    Real :: res
    ! ..
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    ! .. Scalar Arguments ..
    Real :: bound
    Integer :: ii
    ! ..
    ! .. Array Arguments ..
    Real :: qv(:)
    ! ..
    ! .. Local Scalars ..
    Real :: tmp
    Integer :: i
    ! ..
    res = 0.0
    Do i = 1, tp%dim
       tmp = tp%the_data(ii,i) - qv(i)
       res = res + tmp*tmp
       If (res>bound) Then
          res = -1.0
          Return
       End If
    End Do
  End Function bounded_square_distance

  Recursive Subroutine n_search(tp,sr,node)
    ! This is the innermost core routine of the kd-tree search.
    ! it is thus not programmed in quite the clearest style, but is
    ! designed for speed.  -mbk
    ! .. Structure Arguments ..
    Type (tree_node), Pointer :: node
    Type (tree_search_record), Pointer :: sr
    Type (tree_master_record), Pointer :: tp
    ! ..
    ! .. Local Scalars ..
    Real :: dp, sd, sdp, tmp
    Integer :: centeridx, i, ii, j, jmax, k
    Logical :: condition, not_fully_sized
    ! ..
    ! .. Local Arrays ..
    Real, Pointer :: qv(:)
    Integer, Pointer :: ind(:)
    ! ..
    ! ..
    If ( .Not. (ASSOCIATED(node%left)) .And. ( .Not. ASSOCIATED( &
     node%right))) Then
       ! we are on a terminal node
       ind => tp%indexes     ! save for easy access
       qv => sr%qv
       centeridx = sr%centeridx
       ! search through terminal bucket.
       Do i = node%l, node%u
          ii = ind(i)
          condition = .False.
          If (centeridx<0) Then
             condition = .True.
          Else
             condition = (ABS(ii-sr%centeridx)>=sr%correltime)
          End If
          If (condition) Then
             ! sd = bounded_square_distance(tp,qv,ii,sr%bsd)
             ! replace with call to bounded_square distance with inline
             ! code, an
             ! a goto.  Tested to be significantly faster.   SPECIFIC
             ! FOR
             ! the EUCLIDEAN METRIC ONLY! BEWARE!

             sd = 0.0
             Do k = 1, tp%dim
                tmp = tp%the_data(ii,k) - qv(k)
                sd = sd + tmp*tmp
                If (sd>sr%bsd) Then
                   Go To 100
                End If
             End Do
             ! we only consider it if it is better than the 'best' on
             ! the list so far.
             ! if we get here
             ! we know sd is < bsd, and bsd is on the end of the list
             not_fully_sized = (sr%nfound<sr%nbst)
             If (not_fully_sized) Then
                jmax = sr%nfound
                sr%nfound = sr%nfound + 1
             Else
                jmax = sr%nbst - 1
             End If
             ! add it to the list
             Do j = jmax, 1, -1
                If (sd>=sr%dsl(j)) Exit ! we hit insertion location
                sr%il(j+1) = sr%il(j)
                sr%dsl(j+1) = sr%dsl(j)
             End Do
             sr%il(j+1) = ii
             sr%dsl(j+1) = sd
             If ( .Not. not_fully_sized) Then
                sr%bsd = sr%dsl(sr%nbst)
             End If
100          Continue        ! we jumped  here if we had a bad point.
          End If
       End Do
    Else
       ! we are not on a terminal node

       ! Alrighty, this section is essentially the content of the
       ! the Sproul method for searching a Kd-tree, in other words
       ! the second nested "if" statements in the two halves below
       ! and when they get activated.
       dp = sr%qv(node%dnum) - node%val
       If ( .Not. sr%linfinity_metric) Then
          sdp = dp*dp        ! Euclidean
       Else
          sdp = ABS(dp)      ! Linfinity
       End If
       If (dp<0.0) Then
          Call n_search(tp,sr,node%left)
          If (sdp<sr%bsd) Call n_search(tp,sr,node%right)
          ! if the distance projected to the 'wall boundary' is less
          ! than the radius of the ball we must consider, then perform
          ! search on the 'other side' as well.
       Else
          Call n_search(tp,sr,node%right)
          If (sdp<sr%bsd) Call n_search(tp,sr,node%left)
       End If
    End If
  End Subroutine n_search

  Subroutine n_nearest_to_brute_force(tp,qv,n,indexes,distances)
    ! find the 'n' nearest neighbors to 'qv' by exhaustive search.
    ! only use this subroutine for testing, as it is SLOW!  The
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    ! .. Scalar Arguments ..
    Integer, Intent (In) :: n
    ! ..
    ! .. Array Arguments ..
    Real, Target :: distances(n)
    Real, Intent (In) :: qv(:)
    Integer, Target :: indexes(n)
    ! ..
    ! .. Local Scalars ..
    Integer :: i, j, k
    ! ..
    ! .. Local Arrays ..
    Real, Allocatable :: all_distances(:)
    ! ..
    ! ..
    Allocate (all_distances(tp%n))
    Do i = 1, tp%n
       all_distances(i) = bounded_square_distance(tp,qv,i,HUGE(1.0))
    End Do
    ! now find 'n' smallest distances
    Do i = 1, n
       distances(i) = HUGE(1.0)
       indexes(i) = -1
    End Do
    Do i = 1, tp%n
       If (all_distances(i)<distances(n)) Then
          ! insert it somewhere on the list
          Do j = 1, n
             If (all_distances(i)<distances(j)) Exit
          End Do
          ! now we know 'j'
          Do k = n - 1, j, -1
             distances(k+1) = distances(k)
             indexes(k+1) = indexes(k)
          End Do
          distances(j) = all_distances(i)
          indexes(j) = i
       End If
    End Do
    Deallocate (all_distances)
  End Subroutine n_nearest_to_brute_force
End Module kd_tree
!
! MAIN PROGRAM HERE.  This is just an example so you know
! how to use it. 
! 

module time_kdtree
  use kd_tree
contains

  real function time_search(tree,nsearch,mode, nn, r2)
    !
    !  Return CPU time, in seconds, for searching 'nsearch' reference points
    !  using any specific search mode. 
    !

    type(tree_master_record), pointer :: tree
    integer, intent(in)               :: nsearch ! how many reference points
    integer, intent(in)               :: mode    ! what kind of search
    integer, intent(in)               :: nn      ! number of neighbors 
    real, intent(in)                  :: r2      ! radius^2
    !
    real    :: qv(tree%dim), rv               ! query vector, random variate
    integer :: i, random_loc
    real    :: t0, t1
    real, pointer    :: dists(:)
    integer, pointer :: inds(:)
    real(kind(0.0d0)) :: nftotal

    call cpu_time(t0)   ! get initial time.
    nftotal = 0.0
    allocate(dists(nn))
    allocate(inds(nn) )
    do i=1,nsearch

       select case (mode)
       case (1)
          !
          !  Fixed NN search around randomly chosen point
          !
          call random_number(qv)
          call n_nearest_to(tp=tree,qv=qv,n=nn,indexes=inds,distances=dists)
       case (2)
          !
          ! Fixed NN seasrch around randomly chosen point in dataset
          ! with 100 correlation time
          !
          call random_number(rv)
          random_loc = floor(rv*tree%n) + 1
          call n_nearest_to_around_point(tp=tree,idxin=random_loc,&
           correltime=100,n=nn,indexes=inds,distances=dists)
       case default
          write (*,*) 'Search type ', mode, ' not implemented.'
          time_search = -1.0 ! invalid
          return
       end select

    enddo
    call cpu_time(t1)

    time_search = t1-t0
    if (nftotal .gt. 0.0) then
!       write (*,*) 'Average number of neighbors found = ', nftotal / real(nsearch)
    endif
    deallocate(inds,dists)
    return
  end function time_search


  real function searches_per_second(tree,mode,nn,r2) result(res)
    !
    !
    ! return estimated number of searches per second.
    ! Will call "time_search" with increasing numbers of reference points 
    ! until CPU time taken is at least 1 second.
    !

    type(tree_master_record), pointer :: tree
    integer, intent(in)               :: mode    ! what kind of search
    integer, intent(in)               :: nn      ! number of neighbors 
    real, intent(in)                  :: r2      ! radius^2
    !
    integer :: nsearch
    real    :: time_taken

    nsearch = 50  ! start with 50 reference points
    do
       time_taken = time_search(tree,nsearch,mode,nn,r2)
       if (time_taken .lt. 1.0) then
          nsearch = nsearch * 5
          cycle
       else
          res = real(nsearch) / time_taken
          return
       end if
    end do
    return
  end function searches_per_second

end module time_kdtree

Program kd_tree_test
  Use kd_tree  ! this is how you get access to the k-d tree routines
  use time_kdtree

  Integer :: n, d 
  Real, Dimension(:,:), Allocatable :: my_array

  Type(tree_master_record), Pointer :: tree
  ! this is how you declare a tree in your main program

  Integer :: k, j

  Real, Allocatable :: query_vec(:)
  Real,allocatable     :: dists(:), distsb(:) 
  Integer,allocatable   :: inds(:), indsb(:) 
  integer :: nnbrute
  real      :: t0, t1, sps

  integer, parameter  :: nnn = 5
  integer,parameter   :: nnarray(nnn) =(/ 1, 5, 10, 25, 500 /) 
  integer, parameter  :: nr2 = 5

  real      :: r2array(nr2)
  parameter (r2array = (/  1.0e-4,1.0e-3,5.0e-3,1.0e-2,2.0e-2 /) )

  Print *, "Type in N and d"
  Read *, n, d

!  n = 250000
!  d = 3

  Allocate(my_array(n,d))
  allocate(query_vec(d))
  Call Random_number(my_array)  !fills entire array with built-in-randoms

  call cpu_time(t0)
  tree => create_tree(my_array) ! this is how you create a tree. 
  call cpu_time(t1)
  write (*,*) real(n)/real(t1-t0), ' points per second built for regular tree.'

  nnbrute = 50
  allocate(dists(nnbrute),distsb(nnbrute))
  allocate(inds(nnbrute),indsb(nnbrute))

  Do k=1,50
     Call Random_number(query_vec)

     ! find five nearest neighbors to.
     inds = -666
     indsb = -777
     Call n_nearest_to(tp=tree, qv=query_vec, n=nnbrute, indexes=inds, distances=dists)

     call n_nearest_to_brute_force(tree, query_vec, nnbrute, indsb, distsb)
     if (ANY(indsb(1:nnbrute) .ne. inds(1:nnbrute)) .or. &
       any(dists(1:nnbrute) .ne. distsb(1:nnbrute))) then
        write (*,*) 'MISMATCH! @ k=',k
        
        do j=1,nnbrute
           write (*,*) j,'Tree/brute ids=',inds(j),indsb(j)
           write (*,*) j,'tree/brute dists=', dists(j), distsb(j)
!           write (*,*) j,'my dist to tree = ',  &
!            square_distance(d,my_array(:,inds(j)),query_vec)
!           write (*,*) j,'my dist to brute = ',  &
!            square_distance(d,my_array(:,indsb(j)),query_vec)
           
        enddo

        print *, "Tree indexes = ", inds(1:nnbrute)

        print *, "Brute indexes = ", indsb(1:nnbrute)
        print *, "Tree distances = ", dists(1:nnbrute)
        print *, "Brute distances = ", distsb(1:nnbrute)
     endif

  Enddo
  

10 format('R^2 search, r2/d=',G10.2,':',F13.0,A)


20 format(A,' NN=',I7,':',F10.0,' searches/s in ',A)

  do k=1,nnn
     sps =searches_per_second(tree,1,nnarray(k),1.0)
     write (*,20) 'Random pts',nnarray(k),sps,'old regular tree.'
  enddo

  do k=1,nnn
     sps = searches_per_second(tree,2,nnarray(k),1.0)
     write (*,20) 'in-data pts',nnarray(k),sps,'old tree.'
  enddo

  Call destroy_tree(tree)  
  ! this releases memory for the tree BUT NOT THE ARRAY OF DATA YOU PASSED
  ! TO MAKE THE TREE.  

  Deallocate(my_array)
  ! deallocate the memory for the data.

End Program kd_tree_test