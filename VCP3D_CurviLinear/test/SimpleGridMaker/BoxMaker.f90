      program GridMaker
      implicit none

      integer :: ni, nj, nk, i, j, k
      real*8  :: dx, dy, dz, lx, ly, lz, x0, y0, z0

      ! Number of elements in i, j and k-direction
      ni = 18 
      nj = 10
      nk = 10

      lx = 4.5d0  
      ly = 1.0d0
      lz = 2.0d0

      dx = lx/ni
      dy = ly/nj
      dz = lz/nk

      x0 = 0.0d0
      y0 = 0.0d0
      z0 = 0.0d0

      ! Write the grid to fort.10
      i = 1
      write(10) i
      write(10) ni+1, nj+1, nk+1

      write(10) (((i*dx+x0,i=0,ni),j=0,nj),k=0,nk), &
                (((j*dy+y0,i=0,ni),j=0,nj),k=0,nk), &
                (((k*dz+z0,i=0,ni),j=0,nj),k=0,nk)

      end program
