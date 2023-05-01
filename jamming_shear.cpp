#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <algorithm>
#include "BM.h"

#define a1 1.0
#define a2 1.4
#define N1 512
#define N2 512
#define Nl 50
#define phi0 0.70
#define phi_max 0.86
#define dphi_fixed 1.0e-2
#define dphi_a 1.0e-4
#define dphi_g 1.0e-4
#define dt0 1.0e-3
#define dt_max 0.1
#define zeta 1.0
#define dim 2
#define polydispersity 0.0
#define cut 1.4
#define skin 1.0
#define dgamma 1.0e-7

void ini_coord(double (*x)[dim],double L){
  srand((unsigned int)time(NULL));
  for(int i=0;i<N1+N2;i++){
    x[i][0] = (double)rand()/RAND_MAX*L;
    x[i][1] = (double)rand()/RAND_MAX*L;
  }
}

void set_diameter(double *a){
  for(int i=0;i<N1;i++){
    a[i]=a1+polydispersity*gaussian_rand();
  }
  for(int i=0;i<N2;i++){
    a[i+N1]=a2+polydispersity*gaussian_rand();
  }
}

void LE_boundary(double (*x)[dim],double gamma,double L){
  for(int i=0;i<N1+N2;i++){
    if(x[i][0]<x[i][1]*gamma){
      x[i][0]+=L;
    }else if(x[i][0]>x[i][1]*gamma+L){
      x[i][0]-=L;
    }

    if(x[i][1]<0){
      x[i][1]+=L;
      x[i][0]+=L*gamma;
    }else if(x[i][1]>L){
      x[i][1]-=L;
      x[i][0]-=L*gamma;
    }
  }
}

void ini_array(double (*x)[dim]){
  for(int i=0;i<N1+N2;i++)
    for(int j=0;j<dim;j++)
      x[i][j]=0.0;
}

void affine_arithmetic(double (*x)[dim],double *L,double *phi,double dphi){
  int i,j=0;
  // update phi
  *phi+=dphi;
  // transform L and x by using new phi
  *L*=sqrt((*phi-dphi)/(*phi));
  for(i=0;i<N1+N2;i++){
    for(j=0;j<dim;j++){
      x[i][j] *= sqrt((*phi-dphi)/(*phi));
    }
  }
}

void affine_geometric(double (*x)[dim],double *L,double *phi){
  int i,j=0;
  // update phi
  *phi*=(1.0+dphi_g);
  // transform L and x by using old phi
  *L*=sqrt(1.0/(1.0+dphi_g));
  for(i=0;i<N1+N2;i++){
    for(j=0;j<dim;j++){
      x[i][j]*=sqrt(1.0/(1.0+dphi_g));
    }
  }
}

void list_verlet(int (*list)[Nl],double (*x)[dim],double L){
  double dx,dy,dr2;
  int i,j;
  for(i=0;i<N1+N2;i++){
    for(j=0;j<Nl;j++){
      list[i][j]=0;
    }
  }
  
  for(i=0;i<N1+N2;i++){
    for(j=0;j<i;j++){
      dx=x[i][0]-x[j][0];
      dy=x[i][1]-x[j][1];
      dx-=L*floor((dx+0.5*L)/L);
      dy-=L*floor((dy+0.5*L)/L);
      dr2=dx*dx+dy*dy;
      if(dr2<(cut+skin)*(cut+skin)){
        list[i][0]++;
        list[i][(int)list[i][0]]=j;
      }
    }
  }
}

void update(double (*x_update)[dim],double (*x)[dim]){
  for(int i=0;i<N1+N2;i++)
    for(int j=0;j<dim;j++)
      x_update[i][j]=x[i][j];
}

void calc_disp_max(double *disp_max,double (*x)[dim],double (*x_update)[dim],double L){
  double dx,dy;
  double disp;
  for(int i=0;i<N1+N2;i++){
    dx=x[i][0]-x_update[i][0];
    dy=x[i][1]-x_update[i][1];
    dx-=L*floor((dx+0.5*L)/L);
    dy-=L*floor((dy+0.5*L)/L);
    disp = dx*dx+dy*dy;
    if(disp > *disp_max){
      *disp_max =disp;
    }
  }
  //std::cout<<"disp_max is"<<"\t"<<*disp_max<<"\n";
}

void auto_list_update(double *disp_max,double (*x)[dim],double (*x_update)[dim],int (*list)[Nl],double L){
  static int count=0;
  count++;
  calc_disp_max(&(*disp_max),x,x_update,L);
  if(*disp_max > skin*skin*0.25){
    list_verlet(list,x,L);
    update(x_update,x);
    //std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
    *disp_max=0.0;
    count=0;
  }
}

void calc_force(int (*list)[Nl],double (*x)[dim],double (*f)[dim],double *a,double *U,double *p,double *sxy,double gamma,double L){
  double dx,dy,dr,dr2,dUr,aij;
  int i,j;

  *U=0.0;
  *p=0.0;
  *sxy=0.0;
  ini_array(f);
  
  // calculate force
  for(i=0;i<N1+N2;i++){
    for(j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];
      if(dy>0.5*L){
        dy-=L;
        dx-=L*gamma;
      }else if(dy<-0.5*L){
        dy+=L;
        dx+=L*gamma;
      }
      if(dx>0.5*L){
        dx-=L;
      }else if(dx<-0.5*L){
        dx+=L;
      }
      dr2=dx*dx+dy*dy;
      dr=sqrt(dr2);
      aij=0.5*(a[i]+a[list[i][j]]);

      if(dr<aij){
        dUr=-(1.0-dr/aij)/aij;
        f[i][0]-=dUr*dx/dr;
        f[list[i][j]][0]+=dUr*dx/dr;
        f[i][1]-=dUr*dy/dr;
        f[list[i][j]][1]+=dUr*dy/dr;

        *U+=(1.0-dr/aij)*(1.0-dr/aij)/2.0/(double)(N1+N2);  // alpha=2 at Hertzian potential
        *p-=dr*dUr/(2.0*L*L);
        *sxy-=dx*dy*dUr/(2.0*L*L);
      }
    //std::cout<<i<<"\t"<<list[i][j]<<"\t"<<aij-dr<<"\t"<<f[i][0]<<"\t"<<f[i][1]<<std::endl;
    }
  }
}

void apply_shear(int (*list)[Nl],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *sxy,double *disp_max,double gamma,double L){
  for(int i=0;i<N1+N2;i++){
    x[i][0]+=x[i][1]*gamma;
  }
  LE_boundary(x,gamma,L);
  auto_list_update(disp_max,x,x_update,list,L);
  calc_force(list,x,f,a,&(*U),&(*p),&(*sxy),gamma,L);
}

void FIRE(int (*list)[Nl],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *sxy,double *f_tot,double *disp_max,double L){
  double v[N1+N2][dim],P,v_norm,f_norm;
  double dt=dt0,A=0.1;
  int count=0;

  // initialize
  ini_array(v);
  list_verlet(list,x,L);
  calc_force(list,x,f,a,&(*U),&(*p),&(*sxy),0,L);

  for(;;){
    // initialize
    P=0.0;
    *f_tot=0.0;
    // 1. MD integration with Velocity Verlet
    for(int i=0;i<N1+N2;i++){
      for(int j=0;j<dim;j++){
        x[i][j]+=v[i][j]*dt+0.5*f[i][j]*dt*dt;
        v[i][j]+=0.5*f[i][j]*dt;
      }
    }
    calc_force(list,x,f,a,&(*U),&(*p),&(*sxy),0,L);
    for(int i=0;i<N1+N2;i++){
      for(int j=0;j<dim;j++){
        v[i][j]+=0.5*f[i][j]*dt;
      }
    }
    LE_boundary(x,0,L);
    auto_list_update(&(*disp_max),x,x_update,list,L);

    // 3. calculate power P
    for(int i=0;i<N1+N2;i++){
      v_norm=sqrt(v[i][0]*v[i][0]+v[i][1]*v[i][1]);
      f_norm=sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1]);
      *f_tot+=f_norm/(N1+N2);
      P+=f[i][0]*v[i][0]+f[i][1]*v[i][1];
      // 4. velocity modification
      for(int j=0;j<dim;j++){
        v[i][j]=(1.0-A)*v[i][j]+A*f[i][j]*v_norm/(f_norm+DBL_EPSILON);
      }
    }

    // converge criterion
    if(*f_tot<1.0e-12){
      break;
    }

    // 5. evaluate power P
    if(P>=0.0){
      // 6. P>0
      count++;
      if(count>5){
        count=0;
        dt=std::min(1.1*dt,dt_max);
        A*=0.99;
      }
    }else{
      // 7. P<0
      count=0;
      ini_array(v);
      dt*=0.5;
      A=0.1;
    }
  }
}

void output_coord(double (*x)[dim],double *a,double dphi,double gamma){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"coord_FIRE_dphi%.3f_strain%.3e_debug.dat",dphi,gamma);
  file.open(filename);
  for(int i=0;i<N1+N2;i++){
    file <<x[i][0]<<"\t"<<x[i][1]<<"\t"<<a[i]<<std::endl;
  }
  file.close();
}

void output_phiJ(double phi){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"phiJ_FIRE_shear.dat");
  file.open(filename,std::ios::app);
  file << phi << std::endl;
  file.close();
}

void output_pressure(int step,double p,double dphi){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"pressure_FIRE_step%d.dat",step);
  file.open(filename,std::ios::app);
  file << dphi << "\t" << p << std::endl;
  file.close();
}

void output_stress(double sxy0,double sxy,double dphi,double gamma){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"stress_dphi%.3f.dat",dphi);
  file.open(filename,std::ios::app);
  file << gamma << "\t" << abs(sxy-sxy0) << std::endl;
  file.close();
}

void initialize(int (*list)[Nl],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *sxy,double *f_tot,double *phi,double *disp_max,double L){
  *phi=phi0;
  ini_array(x);
  ini_coord(x,L);
  list_verlet(list,x,L);
  FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),L);
}

void make_jamming(int (*list)[Nl],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *sxy,double *f_tot,double *phi,double dphi,double *disp_max,double *L){
  while(*U<1.0e-16){
  //while(phi<=phi_max){
    affine_arithmetic(x,&(*L),&(*phi),dphi);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),*L);
  }
  while(*U>1.0e-12){
    affine_arithmetic(x,&(*L),&(*phi),-dphi);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),*L);
  }
  while(*U>1.0e-14){
    affine_arithmetic(x,&(*L),&(*phi),dphi*1.0e-1);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),*L);
  }
  while(*U>1.0e-16){
    affine_arithmetic(x,&(*L),&(*phi),dphi*1.0e-2);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),*L);
  }

  //output_coord(x,a,0.0,0.0);
  output_phiJ(*phi);
}

int main(){
  double x[N1+N2][dim],x_update[N1+N2][dim],v[N1+N2][dim],f[N1+N2][dim],a[N1+N2];
  double U,p,sxy,sxy0,phi,phi_J,f_tot,gamma;
  double L=sqrt((N1*a1/2.0*a1/2.0+N2*a2/2.0*a2/2.0)*M_PI/phi0),disp_max=0.0;
  char filename[128];
  int list[N1+N2][Nl];
  int k=0,out=1;

  set_diameter(a);
  initialize(list,x,x_update,f,a,&U,&p,&sxy0,&f_tot,&phi,&disp_max,L);

  // calculate phi_J
  make_jamming(list,x,x_update,f,a,&U,&p,&sxy0,&f_tot,&phi,dphi_a,&disp_max,&L);
  phi_J=phi;
  //std::cout<<"phi_J="<<phi_J<<"\t"<<"L="<<L<<std::endl;

  // increase phi
  affine_arithmetic(x,&L,&phi,dphi_fixed);
  list_verlet(list,x,L);
  FIRE(list,x,x_update,f,a,&U,&p,&sxy0,&f_tot,&disp_max,L);
  //output_coord(x,a,dphi_fixed,0.0);
  //std::cout<<"phi="<<phi<<"\t"<<"L="<<L<<std::endl;

  // apply shear
  for(k=0;k<10;k++){  //-7 <= log{10}(phi-phi_J) <= -6
    gamma+=dgamma*1.0;
    apply_shear(list,x,x_update,f,a,&U,&p,&sxy,&disp_max,gamma,L);
    //std::cout<<"gamma="<<gamma<<"\t"<<"|stress|="<<abs(sxy-sxy0)<<std::endl;
    //output_coord(x,a,dphi_fixed,gamma);
    output_stress(sxy0,sxy,dphi_fixed,gamma);

    if(gamma>=1.0){
      gamma=0.0;
      //std::cout<<"reset gamma"<<std::endl;
    }
    FIRE(list,x,x_update,f,a,&U,&p,&sxy,&f_tot,&disp_max,L);
  }
  for(k=0;k<9;k++){   //-6 <= log{10}(phi-phi_J) <= -5
    gamma+=dgamma*1.0e+1;
    apply_shear(list,x,x_update,f,a,&U,&p,&sxy,&disp_max,gamma,L);
    //std::cout<<"gamma="<<gamma<<"\t"<<"|stress|="<<abs(sxy-sxy0)<<std::endl;
    //output_coord(x,a,dphi_fixed,gamma);
    output_stress(sxy0,sxy,dphi_fixed,gamma);

    if(gamma>=1.0){
      gamma=0.0;
      //std::cout<<"reset gamma"<<std::endl;
    }
    FIRE(list,x,x_update,f,a,&U,&p,&sxy,&f_tot,&disp_max,L);
  }
  for(k=0;k<9;k++){   //-5 <= log{10}(phi-phi_J) <= -4
    gamma+=dgamma*1.0e+2;
    apply_shear(list,x,x_update,f,a,&U,&p,&sxy,&disp_max,gamma,L);
    //std::cout<<"gamma="<<gamma<<"\t"<<"|stress|="<<abs(sxy-sxy0)<<std::endl;
    //output_coord(x,a,dphi_fixed,gamma);
    output_stress(sxy0,sxy,dphi_fixed,gamma);

    if(gamma>=1.0){
      gamma=0.0;
      //std::cout<<"reset gamma"<<std::endl;
    }
    FIRE(list,x,x_update,f,a,&U,&p,&sxy,&f_tot,&disp_max,L);
  }
  for(k=0;k<9;k++){   //-4 <= log{10}(phi-phi_J) <= -3
    gamma+=dgamma*1.0e+3;
    apply_shear(list,x,x_update,f,a,&U,&p,&sxy,&disp_max,gamma,L);
    //std::cout<<"gamma="<<gamma<<"\t"<<"|stress|="<<abs(sxy-sxy0)<<std::endl;
    //output_coord(x,a,dphi_fixed,gamma);
    output_stress(sxy0,sxy,dphi_fixed,gamma);

    if(gamma>=1.0){
      gamma=0.0;
      //std::cout<<"reset gamma"<<std::endl;
    }
    FIRE(list,x,x_update,f,a,&U,&p,&sxy,&f_tot,&disp_max,L);
  }
  for(k=0;k<9;k++){   //-3 <= log{10}(phi-phi_J) <= -2
    gamma+=dgamma*1.0e+4;
    apply_shear(list,x,x_update,f,a,&U,&p,&sxy,&disp_max,gamma,L);
    //std::cout<<"gamma="<<gamma<<"\t"<<"|stress|="<<abs(sxy-sxy0)<<std::endl;
    //output_coord(x,a,dphi_fixed,gamma);
    output_stress(sxy0,sxy,dphi_fixed,gamma);

    if(gamma>=1.0){
      gamma=0.0;
      //std::cout<<"reset gamma"<<std::endl;
    }
    FIRE(list,x,x_update,f,a,&U,&p,&sxy,&f_tot,&disp_max,L);
  }
  for(k=0;k<9;k++){   //-2 <= log{10}(phi-phi_J) <= -1
    gamma+=dgamma*1.0e+5;
    apply_shear(list,x,x_update,f,a,&U,&p,&sxy,&disp_max,gamma,L);
    //std::cout<<"gamma="<<gamma<<"\t"<<"|stress|="<<abs(sxy-sxy0)<<std::endl;
    //output_coord(x,a,dphi_fixed,gamma);
    output_stress(sxy0,sxy,dphi_fixed,gamma);

    if(gamma>=1.0){
      gamma=0.0;
      //std::cout<<"reset gamma"<<std::endl;
    }
    FIRE(list,x,x_update,f,a,&U,&p,&sxy,&f_tot,&disp_max,L);
  }
  for(k=0;k<9;k++){   //-1 <= log{10}(phi-phi_J) <= 0
    gamma+=dgamma*1.0e+6;
    apply_shear(list,x,x_update,f,a,&U,&p,&sxy,&disp_max,gamma,L);
    //std::cout<<"gamma="<<gamma<<"\t"<<"|stress|="<<abs(sxy-sxy0)<<std::endl;
    //output_coord(x,a,dphi_fixed,gamma);
    output_stress(sxy0,sxy,dphi_fixed,gamma);
    
    if(gamma>=1.0){
      gamma=0.0;
      std::cout<<"reset gamma"<<std::endl;
    }
    FIRE(list,x,x_update,f,a,&U,&p,&sxy,&f_tot,&disp_max,L);
  }

  return 0;
}
