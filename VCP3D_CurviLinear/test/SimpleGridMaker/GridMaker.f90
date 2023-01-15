      program GridMaker
      implicit none

      integer :: ni, nj, nk, i, j, k
      real*8  :: dx, dy, dz, lx, ly, lz, x0, y0, z0

      ! Number of elements in i, j and k-direction
      ni = 30 
      nj = 15
      nk = 15

      lx = 2.d0  
      ly = 1.d0
      lz = 1.d0

      dx = lx/ni
      dy = ly/nj
      dz = lz/nk

      x0 = -1.0d0
      y0 = -0.5d0
      z0 = -0.5d0

      ! Write the grid to fort.10
      i = 1
      write(10) i
      write(10) ni+1, nj+1, nk+1

      write(10) (((i*dx+x0,i=0,ni),j=0,nj),k=0,nk), &
                (((j*dy+y0,i=0,ni),j=0,nj),k=0,nk), &
                (((k*dz+z0,i=0,ni),j=0,nj),k=0,nk)

      end program
