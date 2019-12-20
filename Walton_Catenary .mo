package Walton_Catenary
  function catThIterator
    input Real l;
    input Real x[:];
    input Real h;
    input Real spm_chain_sub;
    output Real Th[size(x, 1)];
  protected
    constant Real g = Modelica.Constants.g_n;
    Real Thl;
    Real Thlc;
    Real ls[size(x, 1)];
    Real X[size(x, 1)];
  algorithm
    for j in 1:size(x, 1) loop
      if x[j] == 0 then
        Th[j] := 0;
      else
        Thl := 500;
        Thlc := 0;
        while abs(Thlc - Thl) > 0.00001 loop
          Thlc := Thl;
          Th[j] := x[j] * spm_chain_sub * g / Modelica.Math.acosh(1 + spm_chain_sub * g * h / Thlc);
          Thl := Th[j];
        end while;
      end if;
    end for;
  end catThIterator;

  function catXCalculator
    input Real l;
    input Real h;
    input Real spm_chain_sub;
    input Real Th[:];
    input Real x[:];
    output Real X[size(Th, 1)];
  protected
    constant Real g = Modelica.Constants.g_n;
    Real ls[size(Th, 1)];
  algorithm
    for i in 1:size(Th, 1) loop
      if Th[i] == 0 then
        ls[i] := h;
        X[i] := l - ls[i];
      else
        ls[i] := h * (1 + 2 * (Th[i] / (spm_chain_sub * g * h))) ^ 0.5;
        X[i] := l - ls[i] + x[i];
      end if;
    end for;
  end catXCalculator;

  function initCatShape
    input Integer n_seg;
    input Real seg_l;
    input Real spm_wet;
    input Real d;
    input Real X_i;
    input Real h;
    input Real X[:];
    input Real Th[:];
    output Real X_cat[n_seg + 1];
    output Real Z_cat[n_seg + 1];
    output Real lj[n_seg];
  protected
    constant Real g = Modelica.Constants.g_n "Acceleration due to gravity in m/s2";
    parameter Real Th_X_i = linearInterpolatorSV(X, Th, X_i);
    parameter Real a = Th_X_i / (spm_wet * g);
    parameter Real x_i = a * Modelica.Math.acosh(h / a + 1);
    parameter Real X_tdp = X_i - x_i;
    Real s_cat[n_seg + 1];
    Real x_cat[n_seg + 1];
    Real z_cat[n_seg + 1];
    Real s_cat_i;
    Integer count;
    Real X_lnk_whole_floor;
    Real part_s_cat_floor;
  algorithm
    s_cat_i := 0;
    count := 1;
    while s_cat_i <= X_tdp loop
      s_cat[count] := s_cat_i;
      x_cat[count] := s_cat_i;
      X_cat[count] := x_cat[count];
      z_cat[count] := -d;
      Z_cat[count] := z_cat[count];
      s_cat_i := s_cat_i + seg_l;
      X_lnk_whole_floor := X_cat[count];
      count := count + 1;
    end while;
    part_s_cat_floor := X_tdp - X_lnk_whole_floor;
    s_cat_i := 0;
    for i in count:n_seg + 1 loop
      s_cat[i] := seg_l - part_s_cat_floor + s_cat_i;
      x_cat[i] := a * Modelica.Math.asinh(s_cat[i] / a);
      X_cat[i] := X_tdp + x_cat[i];
      z_cat[i] := a * Modelica.Math.cosh(x_cat[i] / a);
      Z_cat[i] := (-d) + z_cat[i] - a;
      s_cat_i := s_cat_i + seg_l;
    end for;
    X_cat[n_seg+1]:=X_i;
    Z_cat[n_seg+1]:=-d+h;
    for i in 1:n_seg loop
      lj[i] := sqrt((X_cat[i + 1] - X_cat[i]) * (X_cat[i + 1] - X_cat[i]) + (Z_cat[i + 1] - Z_cat[i]) * (Z_cat[i + 1] - Z_cat[i]));
    end for;
  end initCatShape;









  function initialTensions
    input Real rho_w;
    input Real d;
    input Real lj_0[n_seg];
    input Integer n_seg;
    input Real spm_dry;
    input Real d_w;
    input Real kj;
    input Real mFu;
    
    input Real xj_c[n_seg + 1];
    input Real zj_c[n_seg + 1];
    input Real dt;
    
    output Real mj_0[n_seg - 1];
    output Real Wj_0[n_seg - 1];
    output Real ejp1_0[n_seg - 1];
    output Real ejm1_0[n_seg - 1];
    
    output Real BetaFu_0;
    output Real EtaFu_0;
  
    output Real Ejm1_c[n_seg];
    output Real Fjm1_c[n_seg];
    output Real Gjm1_c[n_seg];
    output Real Hjm1_c[n_seg];
    
    output Real Ten_c[n_seg];
    output Real xj_n[n_seg + 1];
    output Real zj_n[n_seg + 1];
    
    output Real vxj_n[n_seg + 1];
    output Real vzj_n[n_seg + 1];
    output Real axj_c[n_seg + 1];
    output Real azj_c[n_seg + 1];
    output Real xj_co[n_seg + 1];
    output Real zj_co[n_seg + 1];
  protected
    constant Real g = Modelica.Constants.g_n "Acceleration due to gravity in m/s2";
    constant Real pi = Modelica.Constants.pi "Value of pi";
    parameter Real D = 1.8 * d_w;
    parameter Real csa_m = pi * D * D / 4;
    
    Real cjp1_c[n_seg - 1];
    Real cjm1_c[n_seg - 1];
    Real sjp1_c[n_seg - 1];
    Real sjm1_c[n_seg - 1];
    Real Ij_c[n_seg - 1];
    Real Jj_c[n_seg - 1];
    Real Kj_c[n_seg - 1];
    Real Lj_c[n_seg - 1];
    Real Mj_c[n_seg - 1];
    Real Nj_c[n_seg - 1];
    Real Pj_c[n_seg + 1];
    Real Qj_c[n_seg + 1];
    Real Rj_c[n_seg + 1];
    Real Sj_c[n_seg + 1];
    Real Fu_c[n_seg - 1];
    Real Xj_c[n_seg - 1];
    Real Zj_c[n_seg - 1];
    Real Uj_c[n_seg + 1];
    Real Vj_c[n_seg + 1];
    
    Real Psijm1_c[n_seg];
    Real alpha[n_seg + 1];
    Real beta[n_seg + 1];
    Real den_c;
    Real T05_c;
    Real Ten_tent_c[n_seg];
    Real x_tent_n[n_seg + 1];
    Real z_tent_n[n_seg + 1];
    Real omega_tent_n[n_seg];
    Integer itr_count;
    Real E_tent_n[n_seg];
    Real F_tent_n[n_seg];
    Real G_tent_n[n_seg];
    Real kappa[n_seg + 1];
    Real lamda[n_seg + 1];
    Real den_n;
    Real del_T05_c;
    Real del_Ten_c[n_seg];
    Real del_x_n[n_seg + 1];
    Real del_z_n[n_seg + 1];
    Real lj_n[n_seg];
  algorithm
    BetaFu_0:=8.726547 / (d - mFu * d);
    EtaFu_0:=(-7.25) + BetaFu_0 * d;
    for i in 1:n_seg - 1 loop
      mj_0[i] := 0.5 * (spm_dry * (lj_0[i] + lj_0[i + 1]));
      ejm1_0[i] := rho_w * kj * lj_0[i] * csa_m;
      ejp1_0[i] := rho_w * kj * lj_0[i + 1] * csa_m;
      sjp1_c[i] := (zj_c[i + 2] - zj_c[i + 1]) / lj_0[i + 1];
      cjp1_c[i] := (xj_c[i + 2] - xj_c[i + 1]) / lj_0[i + 1];
      sjm1_c[i] := (zj_c[i + 1] - zj_c[i]) / lj_0[i];
      cjm1_c[i] := (xj_c[i + 1] - xj_c[i]) / lj_0[i];
      Ij_c[i] := mj_0[i] + 0.5 * (ejp1_0[i] * sjp1_c[i] * sjp1_c[i] + ejm1_0[i] * sjm1_c[i] * sjm1_c[i]);
      Jj_c[i] := mj_0[i] + 0.5 * (ejp1_0[i] * cjp1_c[i] * cjp1_c[i] + ejm1_0[i] * cjm1_c[i] * cjm1_c[i]);
      Kj_c[i] := 0.5 * (ejp1_0[i] * sjp1_c[i] * cjp1_c[i] + ejm1_0[i] * sjm1_c[i] * cjm1_c[i]);
      Lj_c[i] := dt * dt * Ij_c[i] / (Ij_c[i] * Jj_c[i] - Kj_c[i] * Kj_c[i]);
      Mj_c[i] := dt * dt * Jj_c[i] / (Ij_c[i] * Jj_c[i] - Kj_c[i] * Kj_c[i]);
      Nj_c[i] := dt * dt * Kj_c[i] / (Ij_c[i] * Jj_c[i] - Kj_c[i] * Kj_c[i]);
      Pj_c[i + 1] := Mj_c[i] * cjm1_c[i] + Nj_c[i] * sjm1_c[i];
      Qj_c[i + 1] := Nj_c[i] * cjm1_c[i] + Lj_c[i] * sjm1_c[i];
      Rj_c[i + 1] := Mj_c[i] * cjp1_c[i] + Nj_c[i] * sjp1_c[i];
      Sj_c[i + 1] := Nj_c[i] * cjp1_c[i] + Lj_c[i] * sjp1_c[i];
      Wj_0[i] := mj_0[i] * g - 0.5 * rho_w * g * (lj_0[i + 1] * csa_m + lj_0[i] * csa_m);
      Xj_c[i] := 0;
      Fu_c[i] := Wj_0[i] * (1 - Modelica.Math.tanh(BetaFu_0 * zj_c[i + 1] + EtaFu_0)) / 2;
      Zj_c[i] := (-Wj_0[i]) + Fu_c[i];
      Uj_c[i + 1] := Mj_c[i] * Xj_c[i] + Nj_c[i] * Zj_c[i];
      Vj_c[i + 1] := Nj_c[i] * Xj_c[i] + Lj_c[i] * Zj_c[i];
    end for;
    for i in 1:n_seg loop
      Ejm1_c[i] := (xj_c[i + 1] - xj_c[i]) * Pj_c[i] + (zj_c[i + 1] - zj_c[i]) * Qj_c[i];
      Fjm1_c[i] := (xj_c[i + 1] - xj_c[i]) * (Pj_c[i + 1] + Rj_c[i]) + (zj_c[i + 1] - zj_c[i]) * (Qj_c[i + 1] + Sj_c[i]);
      Gjm1_c[i] := (xj_c[i + 1] - xj_c[i]) * Rj_c[i + 1] + (zj_c[i + 1] - zj_c[i]) * Sj_c[i + 1];
      Hjm1_c[i] := (xj_c[i + 1] - xj_c[i]) * (Uj_c[i + 1] - Uj_c[i]) + (zj_c[i + 1] - zj_c[i]) * (Vj_c[i + 1] - Vj_c[i]);
    end for;
    Ejm1_c[1] := 0;
    Gjm1_c[n_seg] := 0;
//Psij
    Psijm1_c := Hjm1_c;
//alpha and beta
    alpha[2] := 1;
    beta[2] := 0;
    for i in 3:n_seg + 1 loop
      if abs(Gjm1_c[i - 2]) > 0 then
        alpha[i] := (Fjm1_c[i - 2] * alpha[i - 1] - Ejm1_c[i - 2] * alpha[i - 2]) / Gjm1_c[i - 2];
        beta[i] := (Fjm1_c[i - 2] * beta[i - 1] - Ejm1_c[i - 2] * beta[i - 2] - Psijm1_c[i - 2]) / Gjm1_c[i - 2];
      else
        alpha[i] := 0;
        beta[i] := 0;
      end if;
    end for;
//iterate for tensions
    den_c := Fjm1_c[n_seg] * alpha[n_seg + 1] - Ejm1_c[n_seg] * alpha[n_seg];
    if abs(den_c) < 0.000001 then
      T05_c := 0;
    else
      T05_c := -(Fjm1_c[n_seg] * beta[n_seg + 1] - Ejm1_c[n_seg] * beta[n_seg] - Hjm1_c[n_seg]) / den_c;
    end if;
    Ten_tent_c[1] := T05_c;
    for i in 2:n_seg loop
      Ten_tent_c[i] := alpha[i + 1] * Ten_tent_c[1] + beta[i + 1];
    end for;
    Ten_c := Ten_tent_c;
    omega_tent_n := ones(n_seg);
    itr_count := 0;
    kappa[2] := 1;
    lamda[2] := 0;
    while sum(abs(omega_tent_n)) > 0.01 * n_seg loop
      if itr_count > 1000 then
        terminate("1000 iterations exceeded");
      else
        for i in 1:n_seg + 1 loop
          if i == 1 then
            x_tent_n[i] := xj_c[i];
            z_tent_n[i] := zj_c[i];
          elseif i == n_seg + 1 then
            x_tent_n[i] := xj_c[i] + Uj_c[i];
            z_tent_n[i] := zj_c[i] + Vj_c[i];
          else
            x_tent_n[i] := xj_c[i] + 0.5 * ((-Pj_c[i] * Ten_c[i - 1]) + Rj_c[i] * Ten_c[i] + Uj_c[i]);
            z_tent_n[i] := zj_c[i] + 0.5 * ((-Qj_c[i] * Ten_c[i - 1]) + Sj_c[i] * Ten_c[i] + Vj_c[i]);
          end if;
        end for;
      end if;
      for i in 1:n_seg loop
        omega_tent_n[i] := 0.5 * ((x_tent_n[i + 1] - x_tent_n[i]) * (x_tent_n[i + 1] - x_tent_n[i]) + (z_tent_n[i + 1] - z_tent_n[i]) * (z_tent_n[i + 1] - z_tent_n[i]) - lj_0[i] * lj_0[i]);
        E_tent_n[i] := (x_tent_n[i + 1] - x_tent_n[i]) * Pj_c[i] + (z_tent_n[i + 1] - z_tent_n[i]) * Qj_c[i];
        F_tent_n[i] := (x_tent_n[i + 1] - x_tent_n[i]) * (Pj_c[i + 1] + Rj_c[i]) + (z_tent_n[i + 1] - z_tent_n[i]) * (Qj_c[i + 1] + Sj_c[i]);
        G_tent_n[i] := (x_tent_n[i + 1] - x_tent_n[i]) * Rj_c[i + 1] + (z_tent_n[i + 1] - z_tent_n[i]) * Sj_c[i + 1];
      end for;
      for i in 3:n_seg + 1 loop
        if abs(G_tent_n[i - 2]) > 0 then
          kappa[i] := (F_tent_n[i - 2] * kappa[i - 1] - E_tent_n[i - 2] * kappa[i - 2]) / G_tent_n[i - 2];
          lamda[i] := (F_tent_n[i - 2] * lamda[i - 1] - E_tent_n[i - 2] * lamda[i - 2] - omega_tent_n[i - 2]) / G_tent_n[i - 2];
        else
          kappa[i] := 0;
          lamda[i] := 0;
        end if;
      end for;
      den_n := F_tent_n[n_seg] * kappa[n_seg + 1] - E_tent_n[n_seg] * kappa[n_seg];
      if abs(den_n) < 0.000001 then
        del_T05_c := 0;
      else
        del_T05_c := -(F_tent_n[n_seg] * lamda[n_seg + 1] - E_tent_n[n_seg] * lamda[n_seg] - omega_tent_n[n_seg]) / den_n;
      end if;
      del_Ten_c[1] := del_T05_c;
      for i in 2:n_seg loop
        del_Ten_c[i] := kappa[i + 1] * del_Ten_c[1] + lamda[i + 1];
      end for;
      Ten_c := Ten_c + del_Ten_c;
      itr_count := itr_count + 1;
    end while;
    for i in 2:n_seg loop
      del_x_n[i] := 0.5 * ((-Pj_c[i] * del_Ten_c[i - 1]) + Rj_c[i] * del_Ten_c[i]);
      del_z_n[i] := 0.5 * ((-Qj_c[i] * del_Ten_c[i - 1]) + Sj_c[i] * del_Ten_c[i]);
    end for;
    xj_n := x_tent_n + del_x_n;
    zj_n := z_tent_n + del_z_n;
    for i in 1:n_seg loop
      lj_n[i] := sqrt((xj_n[i + 1] - xj_n[i]) * (xj_n[i + 1] - xj_n[i]) + (zj_n[i + 1] - zj_n[i]) * (zj_n[i + 1] - zj_n[i]));
    end for;
    for i in 1:n_seg + 1 loop
      vxj_n[i] := (xj_n[i] - xj_c[i]) / dt;
      vzj_n[i] := (zj_n[i] - zj_c[i]) / dt;
      axj_c[i] := (xj_n[i] - xj_c[i]) / (dt * dt);
      azj_c[i] := (zj_n[i] - zj_c[i]) / (dt * dt);
    end for;
    xj_co := xj_c;
    zj_co := zj_c;
  end initialTensions;

















  function linearInterpolatorSV
    input Real x[:];
    input Real y[:];
    input Real p;
    output Real q;
  protected
    Integer j;
  algorithm
    if p == x[size(x, 1)] then
      q := y[size(x, 1)];
    else
      j := 1;
      while p >= x[j + 1] loop
        j := j + 1;
      end while;
      q := (y[j + 1] - y[j]) / (x[j + 1] - x[j]) * (p - x[j]) + y[j];
    end if;
  end linearInterpolatorSV;

  function main_WaltonCatenary
    input Real d = 50 "Depth in m";
    input Real rho_w = 1025 "Desnity of water in kg/m3";
    input Real c_max_x=1 "Maximum current velocity";
    input Real seg_l = 10 "Length of mooring segment in m";
    input Integer n_seg = 10 "Number of segments in mooring";
    input Real X_i = 70;
    input Real h = d;
    input Real d_w = 0.022 "Diameter of mooring wire in m";
    input Real cd = 1;
    input Real rho_mat = 7800 "Density of mooring material in kg/m3";
    input Real spm_dry = 10 "Specific mass of mooring in kg/m";
    input Real mFu = 0.999999 "factor to multiply depth where upward thrust on node becomes 5% of full value";
    input Real kj = 1;
    input Real dt = 0.1;
    input Real T_sim=75;
    input Real T_rmp=10;
    input Real T_mov=0;
    input Real ax=4;
    //outputs from initCatShape function
    output Real xj_0[n_seg + 1];
    output Real zj_0[n_seg + 1];
    output Real lj_0[n_seg];
    output Real Ten_0[n_seg];
    //outputs from Tensions functions
    output Real xj_c[n_seg+1];
    output Real zj_c[n_seg+1];
    output Real Ten_c[n_seg];
    output Real T_hist[size(T_s,1)];
    output Real xj_30[n_seg+1];
    output Real zj_30[n_seg+1];
    output Real xj_60[n_seg+1];
    output Real zj_60[n_seg+1];
    output Real z5[size(T_s,1)];
    output Real z4[size(T_s,1)];
    output Real z3[size(T_s,1)];
    output Real z2[size(T_s,1)];
    
  protected
    parameter Real lc = seg_l * n_seg "Length of mooring in m";
    parameter Real spm_wet = spm_dry - spm_dry / rho_mat * rho_w "Submerged weight per meter of chain in kg/m";
    parameter Real X0 = lc - d "Xmin (in m) for iteration of horizontal tension";
    parameter Real Xmax = sqrt(lc ^ 2 - d ^ 2) "Xmax in m";
    parameter Real xmax = Xmax - 10 "Xmax allowable in m";
    parameter Real x[:] = 0:0.1:xmax "Vector of x for iteration of horizontal tension";
    parameter Real Th[size(x, 1)]= catThIterator(lc, x, h, spm_wet);
    parameter Real X[size(x, 1)]= catXCalculator(lc, h, spm_wet, Th, x);
    parameter Real T_s[:]=0:dt:T_sim;
    Integer i;
    
      
      Real t;
      Real c;
      Real Ux;
    //outputs of initialTensions function
      Real mj_0[n_seg - 1];
      Real Wj_0[n_seg - 1];
      Real ejp1_0[n_seg - 1];
      Real ejm1_0[n_seg - 1];
      Real BetaFu_0;
      Real EtaFu_0;
      Real Ejm1_c[n_seg];
      Real Fjm1_c[n_seg];
      Real Gjm1_c[n_seg];
      Real Hjm1_c[n_seg];
      
      Real xj_n[n_seg + 1];
      Real zj_n[n_seg + 1];
      Real lj_n[n_seg];
      Real vxj_n[n_seg + 1];
      Real vzj_n[n_seg + 1];
      Real axj_c[n_seg + 1];
      Real azj_c[n_seg + 1];
   
  algorithm
    
    (xj_0, zj_0, lj_0) := initCatShape(n_seg, seg_l, spm_wet, d, X_i, h, X, Th);
  
    (mj_0,Wj_0,ejp1_0,ejm1_0,BetaFu_0,EtaFu_0,Ejm1_c,Fjm1_c,Gjm1_c,Hjm1_c,Ten_c,xj_n,zj_n,vxj_n,vzj_n,axj_c,azj_c,xj_c,zj_c):=initialTensions(rho_w,d,lj_0,n_seg,spm_dry,d_w,kj,mFu,xj_0,zj_0,dt);
    Ten_0:=Ten_c;
    i:=1;
    T_hist[i]:=Ten_c[n_seg];
    z5[i]:=zj_c[5];
    z4[i]:=zj_c[4];
    z3[i]:=zj_c[3];
    z2[i]:=zj_c[2];
    t:=dt;
    i:=i+1;
  while t<=T_sim+dt loop
    if t<T_rmp then
      c:=c_max_x*t/T_rmp*sin(0.251*t);
    else
      c:=c_max_x*sin(0.251*t);
    end if;
    if t<T_mov then
      Ux:=ax*dt*dt;
    else
      Ux:=0;
    end if;
    (Ten_c,xj_n,zj_n,lj_n,vxj_n,vzj_n,axj_c,azj_c,xj_c,zj_c,Ejm1_c,Fjm1_c,Gjm1_c,Hjm1_c):=tensions(rho_w,c,lj_0,mj_0,Wj_0,ejp1_0,ejm1_0,n_seg,spm_dry,d_w,kj,cd,BetaFu_0,EtaFu_0,xj_n,zj_n,vxj_n,vzj_n,axj_c,azj_c,xj_c,zj_c,Ten_c,Ejm1_c,Fjm1_c,Gjm1_c,Hjm1_c,dt,Ux);
    T_hist[i]:=Ten_c[n_seg];
    z5[i]:=zj_c[5];
    z4[i]:=zj_c[4];
    z3[i]:=zj_c[3];
    z2[i]:=zj_c[2];
    t:=t+dt;
    if i==301 then
    xj_30:=xj_c;
    zj_30:=zj_c;
    end if;
    if i==601 then
    xj_60:=xj_c;
    zj_60:=zj_c;
    end if;
    i:=i+1;
    end while;
end main_WaltonCatenary;










































































































































































































  function tensions
    input Real rho_w;
    input Real c;
    input Real lj_0[n_seg];
    input Real mj_0[n_seg - 1];
    input Real Wj_0[n_seg-1];
    input Real ejp1_0[n_seg - 1];
    input Real ejm1_0[n_seg - 1];
    input Integer n_seg;
    input Real spm_dry;
    input Real d_w;
    input Real kj;
    input Real cd;
    input Real BetaFu_0;
    input Real EtaFu_0;
    input Real xj_c[n_seg + 1];
    input Real zj_c[n_seg + 1];
    input Real vxj_c[n_seg + 1];
    input Real vzj_c[n_seg + 1];
    input Real axj_p[n_seg + 1];
    input Real azj_p[n_seg + 1];
    input Real xj_p[n_seg + 1];
    input Real zj_p[n_seg + 1];
    input Real Ten_p[n_seg];
    input Real Ejm1_p[n_seg];
    input Real Fjm1_p[n_seg];
    input Real Gjm1_p[n_seg];
    input Real Hjm1_p[n_seg];
    input Real dt;
    input Real Ux;
    
    
    output Real Ten_c[n_seg];
    output Real xj_n[n_seg + 1];
    output Real zj_n[n_seg + 1];
    output Real lj_n[n_seg];
    output Real vxj_n[n_seg + 1];
    output Real vzj_n[n_seg + 1];
    output Real axj_c[n_seg + 1];
    output Real azj_c[n_seg + 1];
    output Real xj_co[n_seg + 1];
    output Real zj_co[n_seg + 1];
    output Real Ejm1_c[n_seg];
    output Real Fjm1_c[n_seg];
    output Real Gjm1_c[n_seg];
    output Real Hjm1_c[n_seg];
  protected
    constant Real g = Modelica.Constants.g_n "Acceleration due to gravity in m/s2";
    constant Real pi = Modelica.Constants.pi "Value of pi";
    parameter Real D = 1.8 * d_w;
    parameter Real csa_m = pi * D * D / 4;
    
    Real cjp1_c[n_seg - 1];
    Real cjm1_c[n_seg - 1];
    Real sjp1_c[n_seg - 1];
    Real sjm1_c[n_seg - 1];
    Real Ij_c[n_seg - 1];
    Real Jj_c[n_seg - 1];
    Real Kj_c[n_seg - 1];
    Real Lj_c[n_seg - 1];
    Real Mj_c[n_seg - 1];
    Real Nj_c[n_seg - 1];
    Real Pj_c[n_seg + 1];
    Real Qj_c[n_seg + 1];
    Real Rj_c[n_seg + 1];
    Real Sj_c[n_seg + 1];
    Real qjp1_c[n_seg - 1];
    Real qjm1_c[n_seg - 1];
    Real fjp1_c[n_seg - 1];
    Real fjm1_c[n_seg - 1];
    Real Djp1_c[n_seg - 1];
    Real Djm1_c[n_seg - 1];
    Real Fu_c[n_seg - 1];
    Real Xj_c[n_seg - 1];
    Real Zj_c[n_seg - 1];
    Real Uj_c[n_seg + 1];
    Real Vj_c[n_seg + 1];
    Real Psijm1_c[n_seg];
    Real alpha[n_seg + 1];
    Real beta[n_seg + 1];
    Real den_c;
    Real T05_c;
    Real Ten_tent_c[n_seg];
    Real x_tent_n[n_seg + 1];
    Real z_tent_n[n_seg + 1];
    Real omega_tent_n[n_seg];
    Integer itr_count;
    Real E_tent_n[n_seg];
    Real F_tent_n[n_seg];
    Real G_tent_n[n_seg];
    Real kappa[n_seg + 1];
    Real lamda[n_seg + 1];
    Real den_n;
    Real del_T05_c;
    Real del_Ten_c[n_seg];
    Real del_x_n[n_seg + 1];
    Real del_z_n[n_seg + 1];
  
  algorithm
    for i in 1:n_seg - 1 loop
      sjp1_c[i] := (zj_c[i + 2] - zj_c[i + 1]) / lj_0[i + 1];
      cjp1_c[i] := (xj_c[i + 2] - xj_c[i + 1]) / lj_0[i + 1];
      sjm1_c[i] := (zj_c[i + 1] - zj_c[i]) / lj_0[i];
      cjm1_c[i] := (xj_c[i + 1] - xj_c[i]) / lj_0[i];
      Ij_c[i] := mj_0[i] + 0.5 * (ejp1_0[i] * sjp1_c[i] * sjp1_c[i] + ejm1_0[i] * sjm1_c[i] * sjm1_c[i]);
      Jj_c[i] := mj_0[i] + 0.5 * (ejp1_0[i] * cjp1_c[i] * cjp1_c[i] + ejm1_0[i] * cjm1_c[i] * cjm1_c[i]);
      Kj_c[i] := 0.5 * (ejp1_0[i] * sjp1_c[i] * cjp1_c[i] + ejm1_0[i] * sjm1_c[i] * cjm1_c[i]);
      Lj_c[i] := dt * dt * Ij_c[i] / (Ij_c[i] * Jj_c[i] - Kj_c[i] * Kj_c[i]);
      Mj_c[i] := dt * dt * Jj_c[i] / (Ij_c[i] * Jj_c[i] - Kj_c[i] * Kj_c[i]);
      Nj_c[i] := dt * dt * Kj_c[i] / (Ij_c[i] * Jj_c[i] - Kj_c[i] * Kj_c[i]);
      Pj_c[i + 1] := Mj_c[i] * cjm1_c[i] + Nj_c[i] * sjm1_c[i];
      Qj_c[i + 1] := Nj_c[i] * cjm1_c[i] + Lj_c[i] * sjm1_c[i];
      Rj_c[i + 1] := Mj_c[i] * cjp1_c[i] + Nj_c[i] * sjp1_c[i];
      Sj_c[i + 1] := Nj_c[i] * cjp1_c[i] + Lj_c[i] * sjp1_c[i];
      qjm1_c[i] := (-0.5 * (vxj_c[i + 1] - c + vxj_c[i] - c) * sjm1_c[i]) + 0.5 * (vzj_c[i + 1] + vzj_c[i]) * cjm1_c[i];
      qjp1_c[i] := (-0.5 * (vxj_c[i + 2] - c + vxj_c[i + 1] - c) * sjp1_c[i]) + 0.5 * (vzj_c[i + 2] + vzj_c[i + 1]) * cjp1_c[i];
      fjm1_c[i] := 0.5 * rho_w * cd * D * lj_0[i];
      fjp1_c[i] := 0.5 * rho_w * cd * D * lj_0[i + 1];
      Djm1_c[i] := -fjm1_c[i] * qjm1_c[i] * abs(qjm1_c[i]);
      Djp1_c[i] := -fjp1_c[i] * qjp1_c[i] * abs(qjp1_c[i]);
      
      Fu_c[i] := Wj_0[i] * (1 - Modelica.Math.tanh(BetaFu_0 * zj_c[i + 1] + EtaFu_0)) / 2;
      Xj_c[i] := -0.5 * (Djp1_c[i] * sjp1_c[i] + Djm1_c[i] * sjm1_c[i]);
      Zj_c[i] := 0.5 * (Djp1_c[i] * cjp1_c[i] + Djm1_c[i] * cjm1_c[i]) - Wj_0[i] + Fu_c[i];
      Uj_c[i+1] := Mj_c[i] * Xj_c[i] + Nj_c[i] * Zj_c[i];
      Vj_c[i+1] := Nj_c[i] * Xj_c[i] + Lj_c[i] * Zj_c[i];
    end for;
      Uj_c[n_seg+1]:=Ux;
    for i in 1:n_seg loop
      Ejm1_c[i] := (xj_c[i + 1] - xj_c[i]) * Pj_c[i] + (zj_c[i + 1] - zj_c[i]) * Qj_c[i];
      Fjm1_c[i] := (xj_c[i + 1] - xj_c[i]) * (Pj_c[i + 1] + Rj_c[i]) + (zj_c[i + 1] - zj_c[i]) * (Qj_c[i + 1] + Sj_c[i]);
      Gjm1_c[i] := (xj_c[i + 1] - xj_c[i]) * Rj_c[i + 1] + (zj_c[i + 1] - zj_c[i]) * Sj_c[i + 1];
      Hjm1_c[i] := (xj_c[i + 1] - xj_c[i]) * (Uj_c[i + 1] - Uj_c[i]) + (zj_c[i + 1] - zj_c[i]) * (Vj_c[i + 1] - Vj_c[i]);
    end for;
    Ejm1_c[1] := 0;
    Gjm1_c[n_seg] := 0;
//Psij
    for i in 1:n_seg loop
      if i == 1 then
        Psijm1_c[i] := (-Fjm1_p[i] * Ten_p[i]) + Gjm1_p[i] * Ten_p[i + 1] + Hjm1_p[i] + Hjm1_c[i] + 2 * (xj_c[i + 1] - xj_p[i + 1] - (xj_c[i] - xj_p[i])) * (xj_c[i + 1] - xj_p[i + 1] - (xj_c[i] - xj_p[i])) + 2 * (zj_c[i + 1] - zj_p[i + 1] - (zj_c[i] - zj_p[i])) * (zj_c[i + 1] - zj_p[i + 1] - (zj_c[i] - zj_p[i]));
      elseif i == n_seg then
        Psijm1_c[i] := Ejm1_p[i] * Ten_p[i - 1] - Fjm1_p[i] * Ten_p[i] + Hjm1_p[i] + Hjm1_c[i] + 2 * (xj_c[i + 1] - xj_p[i + 1] - (xj_c[i] - xj_p[i])) * (xj_c[i + 1] - xj_p[i + 1] - (xj_c[i] - xj_p[i])) + 2 * (zj_c[i + 1] - zj_p[i + 1] - (zj_c[i] - zj_p[i])) * (zj_c[i + 1] - zj_p[i + 1] - (zj_c[i] - zj_p[i]));
      else
        Psijm1_c[i] := Ejm1_p[i] * Ten_p[i - 1] - Fjm1_p[i] * Ten_p[i] + Gjm1_p[i] * Ten_p[i + 1] + Hjm1_p[i] + Hjm1_c[i] + 2 * (xj_c[i + 1] - xj_p[i + 1] - (xj_c[i] - xj_p[i])) * (xj_c[i + 1] - xj_p[i + 1] - (xj_c[i] - xj_p[i])) + 2 * (zj_c[i + 1] - zj_p[i + 1] - (zj_c[i] - zj_p[i])) * (zj_c[i + 1] - zj_p[i + 1] - (zj_c[i] - zj_p[i]));
      end if;
    end for;
//alpha and beta
    alpha[2] := 1;
    beta[2] := 0;
    for i in 3:n_seg + 1 loop
      if abs(Gjm1_c[i - 2]) > 0 then
        alpha[i] := (Fjm1_c[i - 2] * alpha[i - 1] - Ejm1_c[i - 2] * alpha[i - 2]) / Gjm1_c[i - 2];
        beta[i] := (Fjm1_c[i - 2] * beta[i - 1] - Ejm1_c[i - 2] * beta[i - 2] - Psijm1_c[i - 2]) / Gjm1_c[i - 2];
      else
        alpha[i] := 0;
        beta[i] := 0;
      end if;
    end for;
//iterate for tensions
    den_c := Fjm1_c[n_seg] * alpha[n_seg + 1] - Ejm1_c[n_seg] * alpha[n_seg];
    if abs(den_c) < 0.000001 then
      T05_c := 0;
    else
      T05_c := -(Fjm1_c[n_seg] * beta[n_seg + 1] - Ejm1_c[n_seg] * beta[n_seg] - Hjm1_c[n_seg]) / den_c;
    end if;
    Ten_tent_c[1] := T05_c;
    for i in 2:n_seg loop
      Ten_tent_c[i] := alpha[i + 1] * Ten_tent_c[1] + beta[i + 1];
    end for;
    Ten_c := Ten_tent_c;
    omega_tent_n := ones(n_seg);
    itr_count := 0;
    kappa[2] := 1;
    lamda[2] := 0;
    while sum(abs(omega_tent_n)) > 0.01 * n_seg loop
      if itr_count > 1000 then
        terminate("1000 iterations exceeded");
      else
        for i in 1:n_seg + 1 loop
          if i == 1 then
            x_tent_n[i] := xj_c[i];
            z_tent_n[i] := zj_c[i];
          elseif i == n_seg + 1 then
            x_tent_n[i] := xj_c[i] + Uj_c[i];
            z_tent_n[i] := zj_c[i] + Vj_c[i];
          else
            x_tent_n[i] := 2*xj_c[i]-xj_p[i] -Pj_c[i] * Ten_c[i - 1] + Rj_c[i] * Ten_c[i] + Uj_c[i];
            z_tent_n[i] := 2*zj_c[i]-zj_p[i] -Qj_c[i] * Ten_c[i - 1] + Sj_c[i] * Ten_c[i] + Vj_c[i];
          end if;
        end for;
      end if;
      for i in 1:n_seg loop
        omega_tent_n[i] := 0.5 * ((x_tent_n[i + 1] - x_tent_n[i]) * (x_tent_n[i + 1] - x_tent_n[i]) + (z_tent_n[i + 1] - z_tent_n[i]) * (z_tent_n[i + 1] - z_tent_n[i]) - lj_0[i] * lj_0[i]);
        E_tent_n[i] := (x_tent_n[i + 1] - x_tent_n[i]) * Pj_c[i] + (z_tent_n[i + 1] - z_tent_n[i]) * Qj_c[i];
        F_tent_n[i] := (x_tent_n[i + 1] - x_tent_n[i]) * (Pj_c[i + 1] + Rj_c[i]) + (z_tent_n[i + 1] - z_tent_n[i]) * (Qj_c[i + 1] + Sj_c[i]);
        G_tent_n[i] := (x_tent_n[i + 1] - x_tent_n[i]) * Rj_c[i + 1] + (z_tent_n[i + 1] - z_tent_n[i]) * Sj_c[i + 1];
      end for;
      for i in 3:n_seg + 1 loop
        if abs(G_tent_n[i - 2]) > 0 then
          kappa[i] := (F_tent_n[i - 2] * kappa[i - 1] - E_tent_n[i - 2] * kappa[i - 2]) / G_tent_n[i - 2];
          lamda[i] := (F_tent_n[i - 2] * lamda[i - 1] - E_tent_n[i - 2] * lamda[i - 2] - omega_tent_n[i - 2]) / G_tent_n[i - 2];
        else
          kappa[i] := 0;
          lamda[i] := 0;
        end if;
      end for;
      den_n := F_tent_n[n_seg] * kappa[n_seg + 1] - E_tent_n[n_seg] * kappa[n_seg];
      if abs(den_n) < 0.000001 then
        del_T05_c := 0;
      else
        del_T05_c := -(F_tent_n[n_seg] * lamda[n_seg + 1] - E_tent_n[n_seg] * lamda[n_seg] - omega_tent_n[n_seg]) / den_n;
      end if;
      del_Ten_c[1] := del_T05_c;
      for i in 2:n_seg loop
        del_Ten_c[i] := kappa[i + 1] * del_Ten_c[1] + lamda[i + 1];
      end for;
      Ten_c := Ten_c + del_Ten_c;
      itr_count := itr_count + 1;
    end while;
    for i in 2:n_seg loop
      del_x_n[i] := -Pj_c[i] * del_Ten_c[i - 1] + Rj_c[i] * del_Ten_c[i];
      del_z_n[i] := -Qj_c[i] * del_Ten_c[i - 1] + Sj_c[i] * del_Ten_c[i];
    end for;
    xj_n := x_tent_n+del_x_n;
    zj_n := z_tent_n+del_z_n;
    for i in 1:n_seg loop
      lj_n[i] := sqrt((xj_n[i + 1] - xj_n[i]) * (xj_n[i + 1] - xj_n[i]) + (zj_n[i + 1] - zj_n[i]) * (zj_n[i + 1] - zj_n[i]));
    end for;
    for i in 1:n_seg + 1 loop
      vxj_n[i] := (xj_n[i] - xj_c[i]) / dt;
      vzj_n[i] := (zj_n[i] - zj_c[i]) / dt;
      axj_c[i] := (xj_n[i] - 2 * xj_c[i] + xj_p[i]) / (dt * dt);
      azj_c[i] := (zj_n[i] - 2 * zj_c[i] + zj_p[i]) / (dt * dt);
    end for;
    xj_co := xj_c;
    zj_co := zj_c;
  end tensions;

  model Fig4_Pseudo
    parameter Integer n_seg = 10;
    Real xj_0[n_seg + 1];
    Real zj_0[n_seg + 1];
    Real lj_0[n_seg];
    Real Ten_0[n_seg];
    Real xj_f[n_seg + 1];
    Real zj_f[n_seg + 1];
    Real Ten_f[n_seg];
   
  equation
    (xj_0, zj_0, lj_0, Ten_0, xj_f, zj_f, Ten_f) = fig4_main_WaltonCatenary();
  end Fig4_Pseudo;























  function fig4_main_WaltonCatenary
  input Real d = 50 "Depth in m";
    input Real rho_w = 1025 "Desnity of water in kg/m3";
    input Real c_max_x = 0 "Maximum current velocity";
    input Real seg_l = 10 "Length of mooring segment in m";
    input Integer n_seg = 10 "Number of segments in mooring";
    input Real X_i = 55;
    input Real h = d;
    input Real d_w = 0.022 "Diameter of mooring wire in m";
    input Real cd = 1;
    input Real rho_mat = 7800 "Density of mooring material in kg/m3";
    input Real spm_dry = 10 "Specific mass of mooring in kg/m";
    input Real mFu = 0.999999 "factor to multiply depth where upward thrust on node becomes 5% of full value";
    input Real kj = 1;
    input Real dt = 0.1;
    input Real T_sim = 0.1;
    input Real T_rmp = 0;
    input Real T_mov = 0;
    input Real ax = 0;
    //outputs from initCatShape function
    output Real xj_0[n_seg + 1];
    output Real zj_0[n_seg + 1];
    output Real lj_0[n_seg];
    output Real Ten_0[n_seg];
    //outputs from Tensions functions
    output Real xj_c[n_seg + 1];
    output Real zj_c[n_seg + 1];
    output Real Ten_c[n_seg];
  protected
    parameter Real lc = seg_l * n_seg "Length of mooring in m";
    parameter Real spm_wet = spm_dry - spm_dry / rho_mat * rho_w "Submerged weight per meter of chain in kg/m";
    parameter Real X0 = lc - d "Xmin (in m) for iteration of horizontal tension";
    parameter Real Xmax = sqrt(lc ^ 2 - d ^ 2) "Xmax in m";
    parameter Real xmax = Xmax - 10 "Xmax allowable in m";
    parameter Real x[:] = 0:0.1:xmax "Vector of x for iteration of horizontal tension";
    parameter Real Th[size(x, 1)] = catThIterator(lc, x, h, spm_wet);
    parameter Real X[size(x, 1)] = catXCalculator(lc, h, spm_wet, Th, x);
    parameter Real T_s[:] = 0:dt:T_sim;
    Real t;
    Real c;
    Real Ux;
    //outputs of initialTensions function
    Real mj_0[n_seg - 1];
    Real Wj_0[n_seg - 1];
    Real ejp1_0[n_seg - 1];
    Real ejm1_0[n_seg - 1];
    Real BetaFu_0;
    Real EtaFu_0;
    Real Ejm1_c[n_seg];
    Real Fjm1_c[n_seg];
    Real Gjm1_c[n_seg];
    Real Hjm1_c[n_seg];
    Real xj_n[n_seg + 1];
    Real zj_n[n_seg + 1];
    Real lj_n[n_seg];
    Real vxj_n[n_seg + 1];
    Real vzj_n[n_seg + 1];
    Real axj_c[n_seg + 1];
    Real azj_c[n_seg + 1];
  algorithm
    (xj_0, zj_0, lj_0) := initCatShape(n_seg, seg_l, spm_wet, d, X_i, h, X, Th);
    (mj_0, Wj_0, ejp1_0, ejm1_0, BetaFu_0, EtaFu_0, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c, Ten_c, xj_n, zj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c) := initialTensions(rho_w, d, lj_0, n_seg, spm_dry, d_w, kj, mFu, xj_0, zj_0, dt);
    Ten_0 := Ten_c;
    t := dt;
    while t <= T_sim + dt loop
      if t < T_rmp then
        c := c_max_x * t / T_rmp;
      else
        c := c_max_x;
      end if;
      if t < T_mov then
        Ux := ax * dt * dt;
      else
        Ux := 0;
      end if;
      (Ten_c, xj_n, zj_n, lj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c) := tensions(rho_w, c, lj_0, mj_0, Wj_0, ejp1_0, ejm1_0, n_seg, spm_dry, d_w, kj, cd, BetaFu_0, EtaFu_0, xj_n, zj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c, Ten_c, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c, dt, Ux);
      t := t + dt;
    end while;
  end fig4_main_WaltonCatenary;

  model Fig5_Pseudo
    parameter Integer n_seg = 10;
    parameter Real T_sim = 50;
    parameter Real dt = 0.1;
    parameter Real T_s[:] = 0:dt:T_sim;
    Real xj_0[n_seg + 1];
    Real zj_0[n_seg + 1];
    Real lj_0[n_seg];
    Real Ten_0[n_seg];
    Real xj_f[n_seg + 1];
    Real zj_f[n_seg + 1];
    Real Ten_f[n_seg];
    Real T_hist[size(T_s, 1)];
  equation
    (xj_0, zj_0, lj_0, Ten_0, xj_f, zj_f, Ten_f, T_hist) = fig5_main_WaltonCatenary();
  end Fig5_Pseudo;









  function fig5_main_WaltonCatenary
  input Real d = 50 "Depth in m";
    input Real rho_w = 1025 "Desnity of water in kg/m3";
    input Real c_max_x = 5 "Maximum current velocity";
    input Real seg_l = 10 "Length of mooring segment in m";
    input Integer n_seg = 10 "Number of segments in mooring";
    input Real X_i = 55;
    input Real h = d;
    input Real d_w = 0.022 "Diameter of mooring wire in m";
    input Real cd = 1;
    input Real rho_mat = 7800 "Density of mooring material in kg/m3";
    input Real spm_dry = 10 "Specific mass of mooring in kg/m";
    input Real mFu = 0.999999 "factor to multiply depth where upward thrust on node becomes 5% of full value";
    input Real kj = 1;
    input Real dt = 0.1;
    input Real T_sim = 50;
    input Real T_rmp = 10;
    input Real T_mov = 0;
    input Real ax = 4;
    //outputs from initCatShape function
    output Real xj_0[n_seg + 1];
    output Real zj_0[n_seg + 1];
    output Real lj_0[n_seg];
    output Real Ten_0[n_seg];
    //outputs from Tensions functions
    output Real xj_c[n_seg + 1];
    output Real zj_c[n_seg + 1];
    output Real Ten_c[n_seg];
    output Real T_hist[size(T_s, 1)];
  protected
    parameter Real lc = seg_l * n_seg "Length of mooring in m";
    parameter Real spm_wet = spm_dry - spm_dry / rho_mat * rho_w "Submerged weight per meter of chain in kg/m";
    parameter Real X0 = lc - d "Xmin (in m) for iteration of horizontal tension";
    parameter Real Xmax = sqrt(lc ^ 2 - d ^ 2) "Xmax in m";
    parameter Real xmax = Xmax - 10 "Xmax allowable in m";
    parameter Real x[:] = 0:0.1:xmax "Vector of x for iteration of horizontal tension";
    parameter Real Th[size(x, 1)] = catThIterator(lc, x, h, spm_wet);
    parameter Real X[size(x, 1)] = catXCalculator(lc, h, spm_wet, Th, x);
    parameter Real T_s[:] = 0:dt:T_sim;
    Integer i;
    Real t;
    Real c;
    Real Ux;
    //outputs of initialTensions function
    Real mj_0[n_seg - 1];
    Real Wj_0[n_seg - 1];
    Real ejp1_0[n_seg - 1];
    Real ejm1_0[n_seg - 1];
    Real BetaFu_0;
    Real EtaFu_0;
    Real Ejm1_c[n_seg];
    Real Fjm1_c[n_seg];
    Real Gjm1_c[n_seg];
    Real Hjm1_c[n_seg];
    Real xj_n[n_seg + 1];
    Real zj_n[n_seg + 1];
    Real lj_n[n_seg];
    Real vxj_n[n_seg + 1];
    Real vzj_n[n_seg + 1];
    Real axj_c[n_seg + 1];
    Real azj_c[n_seg + 1];
  algorithm
    (xj_0, zj_0, lj_0) := initCatShape(n_seg, seg_l, spm_wet, d, X_i, h, X, Th);
    (mj_0, Wj_0, ejp1_0, ejm1_0, BetaFu_0, EtaFu_0, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c, Ten_c, xj_n, zj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c) := initialTensions(rho_w, d, lj_0, n_seg, spm_dry, d_w, kj, mFu, xj_0, zj_0, dt);
    Ten_0 := Ten_c;
    i := 1;
    T_hist[i] := Ten_c[n_seg];
    t := dt;
    i := i + 1;
    while t <= T_sim + dt loop
      if t < T_rmp then
        c := c_max_x * t / T_rmp;
      else
        c := c_max_x;
      end if;
      if t < T_mov then
        Ux := ax * dt * dt;
      else
        Ux := 0;
      end if;
      (Ten_c, xj_n, zj_n, lj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c) := tensions(rho_w, c, lj_0, mj_0, Wj_0, ejp1_0, ejm1_0, n_seg, spm_dry, d_w, kj, cd, BetaFu_0, EtaFu_0, xj_n, zj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c, Ten_c, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c, dt, Ux);
      T_hist[i] := Ten_c[n_seg];
      t := t + dt;
      i := i + 1;
    end while;
  end fig5_main_WaltonCatenary;


  model Fig6_Pseudo
    parameter Integer n_seg = 10;
    parameter Real T_sim = 60;
    parameter Real dt = 0.1;
    parameter Real T_s[:] = 0:dt:T_sim;
    Real xj_0[n_seg + 1];
    Real zj_0[n_seg + 1];
    Real lj_0[n_seg];
    Real Ten_0[n_seg];
    Real xj_f[n_seg + 1];
    Real zj_f[n_seg + 1];
    Real Ten_f[n_seg];
    Real T_hist[size(T_s, 1)];
    Real xj_30[n_seg + 1];
    Real zj_30[n_seg + 1];
  equation
    (xj_0, zj_0, lj_0, Ten_0, xj_f, zj_f, Ten_f, T_hist, xj_30, zj_30) = fig6_main_WaltonCatenary();
  end Fig6_Pseudo;





  function fig6_main_WaltonCatenary
  input Real d = 50 "Depth in m";
    input Real rho_w = 1025 "Desnity of water in kg/m3";
    input Real c_max_x = 0 "Maximum current velocity";
    input Real seg_l = 10 "Length of mooring segment in m";
    input Integer n_seg = 10 "Number of segments in mooring";
    input Real X_i = 55;
    input Real h = d;
    input Real d_w = 0.022 "Diameter of mooring wire in m";
    input Real cd = 1;
    input Real rho_mat = 7800 "Density of mooring material in kg/m3";
    input Real spm_dry = 10 "Specific mass of mooring in kg/m";
    input Real mFu = 0.999999 "factor to multiply depth where upward thrust on node becomes 5% of full value";
    input Real kj = 1;
    input Real dt = 0.1;
    input Real T_sim = 60;
    input Real T_rmp = 0;
    input Real T_mov = 60;
    input Real ax = 4;
    //outputs from initCatShape function
    output Real xj_0[n_seg + 1];
    output Real zj_0[n_seg + 1];
    output Real lj_0[n_seg];
    output Real Ten_0[n_seg];
    //outputs from Tensions functions
    output Real xj_c[n_seg + 1];
    output Real zj_c[n_seg + 1];
    output Real Ten_c[n_seg];
    output Real T_hist[size(T_s, 1)];
    output Real xj_30[n_seg + 1];
    output Real zj_30[n_seg + 1];
  protected
    parameter Real lc = seg_l * n_seg "Length of mooring in m";
    parameter Real spm_wet = spm_dry - spm_dry / rho_mat * rho_w "Submerged weight per meter of chain in kg/m";
    parameter Real X0 = lc - d "Xmin (in m) for iteration of horizontal tension";
    parameter Real Xmax = sqrt(lc ^ 2 - d ^ 2) "Xmax in m";
    parameter Real xmax = Xmax - 10 "Xmax allowable in m";
    parameter Real x[:] = 0:0.1:xmax "Vector of x for iteration of horizontal tension";
    parameter Real Th[size(x, 1)] = catThIterator(lc, x, h, spm_wet);
    parameter Real X[size(x, 1)] = catXCalculator(lc, h, spm_wet, Th, x);
    parameter Real T_s[:] = 0:dt:T_sim;
    Integer i;
    Real t;
    Real c;
    Real Ux;
    //outputs of initialTensions function
    Real mj_0[n_seg - 1];
    Real Wj_0[n_seg - 1];
    Real ejp1_0[n_seg - 1];
    Real ejm1_0[n_seg - 1];
    Real BetaFu_0;
    Real EtaFu_0;
    Real Ejm1_c[n_seg];
    Real Fjm1_c[n_seg];
    Real Gjm1_c[n_seg];
    Real Hjm1_c[n_seg];
    Real xj_n[n_seg + 1];
    Real zj_n[n_seg + 1];
    Real lj_n[n_seg];
    Real vxj_n[n_seg + 1];
    Real vzj_n[n_seg + 1];
    Real axj_c[n_seg + 1];
    Real azj_c[n_seg + 1];
  algorithm
    (xj_0, zj_0, lj_0) := initCatShape(n_seg, seg_l, spm_wet, d, X_i, h, X, Th);
    (mj_0, Wj_0, ejp1_0, ejm1_0, BetaFu_0, EtaFu_0, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c, Ten_c, xj_n, zj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c) := initialTensions(rho_w, d, lj_0, n_seg, spm_dry, d_w, kj, mFu, xj_0, zj_0, dt);
    Ten_0 := Ten_c;
    i := 1;
    T_hist[i] := Ten_c[n_seg];
    t := dt;
    i := i + 1;
    while t <= T_sim + dt loop
      if t < T_rmp then
        c := c_max_x * t / T_rmp;
      else
        c := c_max_x;
      end if;
      if t < T_mov then
        Ux := ax * dt * dt;
      else
        Ux := 0;
      end if;
      (Ten_c, xj_n, zj_n, lj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c) := tensions(rho_w, c, lj_0, mj_0, Wj_0, ejp1_0, ejm1_0, n_seg, spm_dry, d_w, kj, cd, BetaFu_0, EtaFu_0, xj_n, zj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c, Ten_c, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c, dt, Ux);
      T_hist[i] := Ten_c[n_seg];
      t := t + dt;
      if i == 301 then
        xj_30 := xj_c;
        zj_30 := zj_c;
      end if;
      i := i + 1;
    end while;
  end fig6_main_WaltonCatenary;

  model Fig7_Pseudo
    parameter Integer n_seg = 10;
    parameter Real T_sim = 75;
    parameter Real dt = 0.1;
    parameter Real T_s[:] = 0:dt:T_sim;
    Real xj_0[n_seg + 1];
    Real zj_0[n_seg + 1];
    Real lj_0[n_seg];
    Real Ten_0[n_seg];
    Real xj_f[n_seg + 1];
    Real zj_f[n_seg + 1];
    Real Ten_f[n_seg];
    Real T_hist[size(T_s, 1)];
    Real z5[size(T_s, 1)];
    Real z4[size(T_s, 1)];
    Real z3[size(T_s, 1)];
    Real z2[size(T_s, 1)];
    Real xj_30[n_seg + 1];
    Real zj_30[n_seg + 1];
    Real xj_60[n_seg + 1];
    Real zj_60[n_seg + 1];
  equation
    (xj_0, zj_0, lj_0, Ten_0, xj_f, zj_f, Ten_f, T_hist, xj_30, zj_30, xj_60, zj_60, z5, z4, z3, z2) = main_WaltonCatenary();
  end Fig7_Pseudo;






  function fig7_main_WaltonCatenary
  input Real d = 50 "Depth in m";
    input Real rho_w = 1025 "Desnity of water in kg/m3";
    input Real c_max_x = 1 "Maximum current velocity";
    input Real seg_l = 10 "Length of mooring segment in m";
    input Integer n_seg = 10 "Number of segments in mooring";
    input Real X_i = 70;
    input Real h = d;
    input Real d_w = 0.022 "Diameter of mooring wire in m";
    input Real cd = 1;
    input Real rho_mat = 7800 "Density of mooring material in kg/m3";
    input Real spm_dry = 10 "Specific mass of mooring in kg/m";
    input Real mFu = 0.999999 "factor to multiply depth where upward thrust on node becomes 5% of full value";
    input Real kj = 1;
    input Real dt = 0.1;
    input Real T_sim = 75;
    input Real T_rmp = 10;
    input Real T_mov = 0;
    input Real ax = 4;
    //outputs from initCatShape function
    output Real xj_0[n_seg + 1];
    output Real zj_0[n_seg + 1];
    output Real lj_0[n_seg];
    output Real Ten_0[n_seg];
    //outputs from Tensions functions
    output Real xj_c[n_seg + 1];
    output Real zj_c[n_seg + 1];
    output Real Ten_c[n_seg];
    output Real T_hist[size(T_s, 1)];
    output Real xj_30[n_seg + 1];
    output Real zj_30[n_seg + 1];
    output Real xj_60[n_seg + 1];
    output Real zj_60[n_seg + 1];
    output Real z5[size(T_s, 1)];
    output Real z4[size(T_s, 1)];
    output Real z3[size(T_s, 1)];
    output Real z2[size(T_s, 1)];
  protected
    parameter Real lc = seg_l * n_seg "Length of mooring in m";
    parameter Real spm_wet = spm_dry - spm_dry / rho_mat * rho_w "Submerged weight per meter of chain in kg/m";
    parameter Real X0 = lc - d "Xmin (in m) for iteration of horizontal tension";
    parameter Real Xmax = sqrt(lc ^ 2 - d ^ 2) "Xmax in m";
    parameter Real xmax = Xmax - 10 "Xmax allowable in m";
    parameter Real x[:] = 0:0.1:xmax "Vector of x for iteration of horizontal tension";
    parameter Real Th[size(x, 1)] = catThIterator(lc, x, h, spm_wet);
    parameter Real X[size(x, 1)] = catXCalculator(lc, h, spm_wet, Th, x);
    parameter Real T_s[:] = 0:dt:T_sim;
    Integer i;
    Real t;
    Real c;
    Real Ux;
    //outputs of initialTensions function
    Real mj_0[n_seg - 1];
    Real Wj_0[n_seg - 1];
    Real ejp1_0[n_seg - 1];
    Real ejm1_0[n_seg - 1];
    Real BetaFu_0;
    Real EtaFu_0;
    Real Ejm1_c[n_seg];
    Real Fjm1_c[n_seg];
    Real Gjm1_c[n_seg];
    Real Hjm1_c[n_seg];
    Real xj_n[n_seg + 1];
    Real zj_n[n_seg + 1];
    Real lj_n[n_seg];
    Real vxj_n[n_seg + 1];
    Real vzj_n[n_seg + 1];
    Real axj_c[n_seg + 1];
    Real azj_c[n_seg + 1];
  algorithm
    (xj_0, zj_0, lj_0) := initCatShape(n_seg, seg_l, spm_wet, d, X_i, h, X, Th);
    (mj_0, Wj_0, ejp1_0, ejm1_0, BetaFu_0, EtaFu_0, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c, Ten_c, xj_n, zj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c) := initialTensions(rho_w, d, lj_0, n_seg, spm_dry, d_w, kj, mFu, xj_0, zj_0, dt);
    Ten_0 := Ten_c;
    i := 1;
    T_hist[i] := Ten_c[n_seg];
    z5[i] := zj_c[5];
    z4[i] := zj_c[4];
    z3[i] := zj_c[3];
    z2[i] := zj_c[2];
    t := dt;
    i := i + 1;
    while t <= T_sim + dt loop
      if t < T_rmp then
        c := c_max_x * t / T_rmp * sin(0.251 * t);
      else
        c := c_max_x * sin(0.251 * t);
      end if;
      if t < T_mov then
        Ux := ax * dt * dt;
      else
        Ux := 0;
      end if;
      (Ten_c, xj_n, zj_n, lj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c) := tensions(rho_w, c, lj_0, mj_0, Wj_0, ejp1_0, ejm1_0, n_seg, spm_dry, d_w, kj, cd, BetaFu_0, EtaFu_0, xj_n, zj_n, vxj_n, vzj_n, axj_c, azj_c, xj_c, zj_c, Ten_c, Ejm1_c, Fjm1_c, Gjm1_c, Hjm1_c, dt, Ux);
      T_hist[i] := Ten_c[n_seg];
      z5[i] := zj_c[5];
      z4[i] := zj_c[4];
      z3[i] := zj_c[3];
      z2[i] := zj_c[2];
      t := t + dt;
      if i == 301 then
        xj_30 := xj_c;
        zj_30 := zj_c;
      end if;
      if i == 601 then
        xj_60 := xj_c;
        zj_60 := zj_c;
      end if;
      i := i + 1;
    end while;
  end fig7_main_WaltonCatenary;








































end Walton_Catenary;
