module high_order
  use constants, only: i2,i4,i8,r4,r8

  use datastruct
    use smeshpack  

  contains
  
    subroutine finite_volume_high_order(mesh)

      implicit none      
      
      integer     :: i,j
      real(r8)    :: r
      integer(i4) :: n
      real(r8)    :: pt(1:3)
      integer(i4),allocatable   :: listv(:)
      real(r8),allocatable      :: listd(:)
      character(len=100)        :: fmt
      
      type(grid_structure)      :: mesh
    
      r = 0
      n = 8

      do i = 1,mesh%nv
        write(fmt,*)n
        pt(1) = mesh%v(i)%p(1)
        pt(2) = mesh%v(i)%p(2)
        pt(3) = mesh%v(i)%p(3)  
        write(*,"(3X,'i',2X,'mesh%v(i)%p(1)',2x,'mesh%v(i)%p(2)',2x,'mesh%v(i)%p(3)',4X,'n na entrada')")
        write(*,"(I4,3(5XF9.5),10X,I4)")i,mesh%v(i)%p(1),mesh%v(i)%p(2),mesh%v(i)%p(3),n
        call getnearnodes(pt, mesh, r, listv, listd, n)
        write(*,*)
        write(*,"(2X,'Vizinhos',20X,'n na saida')")
        write(*,"("//adjustl(fmt)//"I4,10X,I4)")(listv(j),j=1,n),n
        pause       
      end do
      return
    end subroutine finite_volume_high_order
end module high_order
