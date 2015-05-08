function [DATA, ns, nt] =  Model(GW,Nens)   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %           SLAM: Channel Version.
    %
    %   Single Level Atmospheric Model
    %   based on the Shallow Water Equations.
    %
    %   u is the zonal velocity component 
    %   v is the meridional velocity component 
    %   z is the geopotential height.
    %
    %   Shallow Water Equations:
    %         du/dt -  fv +  g dz/dx = 0
    %         dv/dt +  fu +  g dz/dy = 0
    %         dz/dt + z(du/dx+dv/dy) = 0
    %
    %   A-grid: (u,v,z) given at all points.
    %   (This causes noise, which must be controlled)
    %
    %   Timestep: Leapfrog (first step forward)
    %
    %   Initial Conditions:
    %                       z0 = zRef + zBar + zWave
    %                       zRef = 5500 metres
    %                       zBar for mean zonal flow
    %                       zWave for perturbation
    %
    %         Initial winds (u0,v0) are geostrophic.
    % 
    %   Boundary Conditions:  Periodic in x
    %                         Walls at N and S
    %                         (v = 0, dzdx = 0)
    % 
    %   Damping:
    %            Second-order diffusion, equivalent to
    %            additional terms in the equations, e.g.
    %            du/dt -  fv +  g dz/dx  =  nu del^2 u 
    %
    %   Fourier Filtering: 
    %           Even with second-order diffusion, there is
    %           a severe problem with two-grid-interval noise.
    %           This is rectified by Fourier-clipping the
    %           smallest scales. The extent of the clipping
    %           is determined by the parameter Fband.
    %     
    %   GW is the wind parameter
    %   
    %   Nens is the number of outputs
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Define basic parameters.
    a = (4e7)/(2*pi);            %  Radius of the Earth
    Omega = 2*pi / ( 60*60*24);  %  Rotation rate of the Earth
    g = 9.8066;                  %  Acceleration of gravity
    nu = 1.25*10^6;              %  Damping coefficient

    Filter = 1;   %  (1) Fourier Filtering (0) No Filtering.
    Fband = 1;    %  Parameter for Filtering band-width.
     
    % Parameters specifying the domain (4:1 aspect ratio).
    global Lx Ly Nx Ny;
    Lx = 32*10^6;   %  Channel length in x-direction
    Ly = 8*10^6;    %  Channel width in y-direction
    Nx = 41;        %  Number of points in x-direction
    Ny = 11;        %  Number of points in y-direction
 
    y = linspace(0,Ly,Ny);   %  Independent variable y
    x = linspace(0,Lx,Nx);   %  Independent variable x
    [X, Y] = meshgrid(x,y);  %  2D arrays, mainly for plotting
    X = X'; Y = Y';          %  MatLab is strange!
    
    dx = Lx/(Nx-1);    %  Grid step in x-direction
    dy = Ly/(Ny-1);    %  Grid step in y-direction
    
    % Set the forcast length and time step
    Tdays = 3;              % Forecast length in days
    Thours = 96*Tdays;      % Forecast length in hours
    Tsecs = Thours*3600;    % Forecast length in seconds
    dt = 1800;              %  Timestep in seconds
    Nsteps = Tsecs / dt;   %  Total number of timesteps 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coriolis parameter (beta-plane version).
    f0 = 2*Omega*sin(pi/4);        %  Mean value of f 
    beta = 2*Omega*cos(pi/4)/a;    %  Constant beta-parameter
    y0 = Ly/2;                     %  Central value of y
    f = f0 + beta*(Y-y0);          %  Coriolis parameter.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Define initial conditions
    
    mx = 1;      %  Wave number in x-direction
    my = 1;      %  Wave number in y-direction

    zRef = 5500;  %  Constant reference geopotential height

    % Initial geopotential height: cos-squared in y-direction.
    y0 = Ly/2;
    zAmp = 500;
    zWave = zAmp*sin(2*pi*mx*X/Lx).*cos(my*pi*(Y-y0)/Ly);
     
    %  Add a mean flow (if uBar non-zero).
    uBar = 100;   uBar = 0;

    %  Mean height from geopotential balance with uBar.
    zBar = -(uBar/g)*(f0*(Y-y0)+0.5*beta*(Y-y0).^2);
    
    %  Construct the initial geopotential height
    z0 = zRef + zWave + zBar;

    ns = numel(z0);
    nt  = Nsteps+1;
    DATA = zeros(ns*nt,Nens);
    
    
    for CI = 1:Nens
    %  Store the geopotential height.
    
    % Perturb the geopotential height component of the initial condition
    perturb=randn(size(z0));
    perturb_coeff = 500;

    z = z0.*(1+perturb/perturb_coeff);
    
    %  Plot initial geopotential height z
    SURF = 0;
    if (SURF==1 )
      surf(X,Y,z); shading('interp');
      axis([0,Lx,0,Ly]);  %  ,5000,6000]);
      view(0,90);
      xlabel('X'); ylabel('Y'); zlabel('Z')
      title('Initial Geopotential Height')
      drawnow
    else
      lines = 4800:40:6200;
      %contourf(X,Y,z,lines)
      %drawnow
    end
    
    %  Calculate the x-gradient of z
    for jj=1:Ny
      for ii=2:Nx-1
        dzdx(ii,jj) = (z(ii+1,jj)-z(ii-1,jj))/(2*dx);
      end
      dzdx( 1,jj) = (z(2,jj)-z(Nx-1,jj))/(2*dx);
      dzdx(Nx,jj) = dzdx(1,jj);
    end
    % Set N and S boundary x-gradients to zero explicitly.
    dzdx(:,1) = 0; dzdx(:,Ny) = 0; 
    
    %  Calculate the y-gradient of z
    for ii=1:Nx
      for jj=2:Ny-1
        dzdy(ii,jj) = (z(ii,jj+1)-z(ii,jj-1))/(2*dy);
      end
      dzdy(ii,1)  = (z(ii,2) -z(ii,1)   )/(dy);
      dzdy(ii,Ny) = (z(ii,Ny)-z(ii,Ny-1))/(dy);
    end
    
    %  Define geostrophic initial winds.
    u0 = -g*dzdy./f;
    v0 = +g*dzdx./f;

    %  Store the initial winds
    u = u0;
    v = v0;

    %  Modify ICs to accentuate Gravity Waves.
    if(GW==1)                    end    %  No action.
    if(GW==2)  u = 0*u; v = 0*v; end    %  Zero winds.
    if(GW==3)  u =  -u; v =  -v; end    %  Reverse winds.

    %  Store the fields for start-up of leapfrog scheme.
    uOld = u;
    vOld = v;
    zOld = z;

    % Compute the Rossby wave-speed
    kx = (2*pi)/(Lx/mx); ky = (2*pi)/(Ly/my); 
    cR = uBar - beta/(kx^2+ky^2);
    Period = Lx/(mx*cR)/(24*3600);

    %  Define the sizes of the second derivative fields.
    d2zdx2 = zeros(Nx,Ny);
    d2udx2 = zeros(Nx,Ny);
    d2vdx2 = zeros(Nx,Ny);
    d2zdy2 = zeros(Nx,Ny);
    d2udy2 = zeros(Nx,Ny);
    d2vdy2 = zeros(Nx,Ny);

    %  Define the sizes of the Laplacians
    lapz = zeros(Nx,Ny);
    lapu = zeros(Nx,Ny);
    lapv = zeros(Nx,Ny);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    Main time-stepping loop.
    
    Z = zeros(numel(z),Nsteps+1);
    
    for nn=1:Nsteps
    
        % Store data
        Z(:,nn) = reshape(z,numel(z),1);
        
    %  Calculate the x-gradients
    for jj=1:Ny
        for ii=2:Nx-1
            dzdx(ii,jj) = (z(ii+1,jj)-z(ii-1,jj))/(2*dx);
            dudx(ii,jj) = (u(ii+1,jj)-u(ii-1,jj))/(2*dx);
            dvdx(ii,jj) = (v(ii+1,jj)-v(ii-1,jj))/(2*dx);
        end
        dzdx( 1,jj) = (z(2,jj)-z(Nx-1,jj))/(2*dx);
        dudx( 1,jj) = (u(2,jj)-u(Nx-1,jj))/(2*dx);
        dvdx( 1,jj) = (v(2,jj)-v(Nx-1,jj))/(2*dx);
        dzdx(Nx,jj) = dzdx(1,jj);
        dudx(Nx,jj) = dudx(1,jj);
        dvdx(Nx,jj) = dvdx(1,jj);
    end
    % Set N and S boundary values to zero explicitly.
    dzdx(:,1) = 0; dzdx(:,Ny) = 0; 
    dvdx(:,1) = 0; dvdx(:,Ny) = 0; 
         
    %  Calculate the y-gradient
    for ii=1:Nx
        for jj=2:Ny-1
          dzdy(ii,jj) = (z(ii,jj+1)-z(ii,jj-1))/(2*dy);
          dudy(ii,jj) = (u(ii,jj+1)-u(ii,jj-1))/(2*dy);
          dvdy(ii,jj) = (v(ii,jj+1)-v(ii,jj-1))/(2*dy);
        end
        dzdy(ii, 1) = ( z(ii, 2) - z0(ii,   1))/(dy);
        dudy(ii, 1) = ( u(ii, 2) -  u(ii,   1))/(dy);
        dvdy(ii, 1) = ( v(ii, 2) -       0    )/(dy);
        dzdy(ii,Ny) = (z0(ii,Ny) -  z(ii,Ny-1))/(dy);
        dudy(ii,Ny) = ( u(ii,Ny) -  u(ii,Ny-1))/(dy);
        dvdy(ii,Ny) = (    0     -  v(ii,Ny-1))/(dy);
    end
    
    % Compute the divergence.
    div = dudx + dvdy;
    
    %%  Calculate the Right-hand sides of the equations.
    Ru = -(u.*dudx + v.*dudy  - f.*v + g*dzdx);
    Rv = -(u.*dvdx + v.*dvdy  + f.*u + g*dzdy);
    Rz = -(u.*dzdx + v.*dzdy  +    z.*div    );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Add the damping term to the Right-hand sides
    Ru = Ru + nu*lapu;
    Rv = Rv + nu*lapv;
    Rz = Rz + nu*lapz;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    STEP FORWARD IN TIME

    % Set the time step
    tstep = 2*dt;
    if(nn==1)
        tstep = dt; % Forward step to begin.
    end
    
    % Step the equations forward:
    % Compute the values at (nn+1)*dt.
    uNew = uOld + tstep*Ru;
    vNew = vOld + tstep*Rv;
    zNew = zOld + tstep*Rz;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Compute damping using the fields at (n-1)dt.
    %
    %   Calculate the Laplacians of z, u and v
    for ii=2:Nx-1
    for jj=2:Ny-1
      d2zdx2(ii,jj) = (z(ii+1,jj)+z(ii-1,jj)-2*z(ii,jj))/(dx^2);
      d2udx2(ii,jj) = (u(ii+1,jj)+u(ii-1,jj)-2*u(ii,jj))/(dx^2);
      d2vdx2(ii,jj) = (v(ii+1,jj)+v(ii-1,jj)-2*v(ii,jj))/(dx^2);

      d2zdy2(ii,jj) = (z(ii,jj+1)+z(ii,jj-1)-2*z(ii,jj))/(dy^2);
      d2udy2(ii,jj) = (u(ii,jj+1)+u(ii,jj-1)-2*u(ii,jj))/(dy^2);
      d2vdy2(ii,jj) = (v(ii,jj+1)+v(ii,jj-1)-2*v(ii,jj))/(dy^2);
    end
    end
    for jj=2:Ny-1
      d2zdx2(1,jj) = (z(2,jj)+z(Nx-1,jj)-2*z(1,jj))/(dx^2);
      d2udx2(1,jj) = (u(2,jj)+u(Nx-1,jj)-2*u(1,jj))/(dx^2);
      d2vdx2(1,jj) = (v(2,jj)+v(Nx-1,jj)-2*v(1,jj))/(dx^2);

      d2zdy2(1,jj) = (z(1,jj+1)+z(1,jj-1)-2*z(1,jj))/(dy^2);
      d2udy2(1,jj) = (u(1,jj+1)+u(1,jj-1)-2*u(1,jj))/(dy^2);
      d2vdy2(1,jj) = (v(1,jj+1)+v(1,jj-1)-2*v(1,jj))/(dy^2);

      d2zdx2(Nx,jj) = d2zdx2(1,jj); d2zdy2(Nx,jj) = d2zdy2(1,jj);
      d2udx2(Nx,jj) = d2udx2(1,jj); d2udy2(Nx,jj) = d2udy2(1,jj);
      d2vdx2(Nx,jj) = d2vdx2(1,jj); d2vdy2(Nx,jj) = d2vdy2(1,jj);
    end

    %  Now get the Laplacians (for use next time round)
    lapz = d2zdx2 + d2zdy2;
    lapu = d2udx2 + d2udy2;
    lapv = d2vdx2 + d2vdy2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Move the nn-fields to (nn-1)-slots
    uOld = u;
    vOld = v;
    zOld = z;
    
    %  Move the (nn+1)-fields to (nn)-slots
    u = uNew;
    v = vNew;
    z = zNew;
    
    %  Store some  central values
    uMid(nn) = uNew((Nx+1)/2,(Ny+1)/2);
    vMid(nn) = vNew((Nx+1)/2,(Ny+1)/2);
    zMid(nn) = zNew((Nx+1)/2,(Ny+1)/2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % Force geostrophic u at N and S boundaries.
    % u(:, 1) = -( (Fi(:, 2)-Fi(:,   1))/dy )./f(:, 1);
    % u(:,Ny) = -( (Fi(:,Ny)-Fi(:,Ny-1))/dy )./f(:,Ny);

    % Force zero divergence at the N and S boundaries.

    % First, apply consistency condition, that the integral
    % of v along the second row up or down vanishes.
    vSodd =  mean(v(1:2:Nx,2));
    vSevn =  mean(v(2:2:Nx,2));
    v(1:2:Nx,2) = v(1:2:Nx,2) - vSodd;
    v(2:2:Nx,2) = v(2:2:Nx,2) - vSevn;
    vNodd =  mean(v(1:2:Nx,Ny-1));
    vNevn =  mean(v(2:2:Nx,Ny-1));
    v(1:2:Nx,Ny-1) = v(1:2:Nx,Ny-1) - vNodd;
    v(2:2:Nx,Ny-1) = v(2:2:Nx,Ny-1) - vNevn;

    % Now integrate along boundaries.
    % u(2, 1) = u(Nx-1, 1) - (2*dx/dy)*v(1,   2);
    % u(2,Ny) = u(Nx-1,Ny) + (2*dx/dy)*v(1,Ny-1);
    for ii=2:Nx-1
      u(ii+1, 1) = u(ii-1, 1) - (2*dx/dy)*v(ii,   2);
      u(ii+1,Ny) = u(ii-1,Ny) + (2*dx/dy)*v(ii,Ny-1);
    end
    % u(1, 1) = u(Nx, 1); 
    % u(1,Ny) = u(Nx,Ny); 

    % Adjust u to maintain average value on boundaries.
    uS = mean(u0(:, 1));
    u(:, 1) = u(:, 1)+uS;
    uN = mean(u0(:,Ny));
    u(:,Ny) = u(:,Ny)+uN;
    
    % Re-set S and N boundary values for z and v.
    z(:,1) = z0(:,1); z(:,Ny) = z0(:,Ny);
    v(:,1) =     0  ;  v(:,Ny) =    0   ;

    % Re-set S and N boundary values for u.
     u(:,1) = u0(:,1) ;  u(:,Ny) = u0(:,Ny) ;
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the z field each timestep.
  
    SURF = 0 ; 
    if (SURF==1 )
      figure(1)
      subplot(1,2,1)
       surf(X,Y,z); shading('interp');
       axis([0,Lx,0,Ly ,5000,6000]);
      view(-9,30);  
     subplot(1,2,2)
      R=5500+(z-5500)*4e0; H=Y; Theta=X/max(max(X))*2*pi; [X1,Y1,Z1] = pol2cart(Theta,R,H); 
        surf(X1,Z1,Y1);shading('faceted')
      view(140,20);  
axis([-7.5e3, 7.5e3, 0, 8e6 -7e3, 7e3 ])
      xlabel('X'); ylabel('Y'); zlabel('Z');
      title('Geopotential Height')
%     colorbar
      drawnow
    else
      %contourf(X,Y,z,lines)
      %drawnow
    end

   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        Fourier Filtering in space.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(Filter==1) 

      %  compute the s-to-N trend and remove it.
      for nx=1:Nx
        zTrend(nx,:) = z(nx,1) + (y/Ly).*(z(nx,Ny)-z(nx,1));
        uTrend(nx,:) = u(nx,1) + (y/Ly).*(u(nx,Ny)-u(nx,1));
        vTrend(nx,:) = v(nx,1) + (y/Ly).*(v(nx,Ny)-v(nx,1));
      end
      zDash = z - zTrend;
      uDash = u - uTrend;
      vDash = v - vTrend;
      
      %  Define filtering masks in x and y directions.
      NxHalf = (Nx+1)/2;
      NyHalf = (Ny+1)/2;
      xMask = ones(Nx,Ny);
      yMask = ones(Nx,Ny);

      xMask(NxHalf-Fband:NxHalf+Fband,:) = 0;
      yMask(:,NyHalf-Fband:NyHalf+Fband) = 0;
      Mask = xMask.*yMask;
      
      zHat = fft2(zDash); zHatFilt = zHat.*Mask;
      uHat = fft2(uDash); uHatFilt = uHat.*Mask;
      vHat = fft2(vDash); vHatFilt = vHat.*Mask;
      
      zDashFilt = real(ifft2(zHatFilt));
      uDashFilt = real(ifft2(uHatFilt));
      vDashFilt = real(ifft2(vHatFilt));
      
      for nx=1:Nx
        z(nx,:) = zDashFilt(nx,:)+zTrend(nx,:);
        u(nx,:) = uDashFilt(nx,:)+uTrend(nx,:);
        v(nx,:) = vDashFilt(nx,:)+vTrend(nx,:);
      end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

              
    end % End of time-stepping loop.

        % Store data
        Z(:,Nsteps+1) = reshape(z,numel(z),1);
        DATA(:,CI) = reshape(Z,numel(Z),1);
    end
    
   

        
        
return
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
