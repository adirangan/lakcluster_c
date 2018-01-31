void hsv2rgb(double h,double s,double v,double *r,double *g,double *b)
{
  int i;
  double f=0,p=0,q=0,t=0;
  if (s==0){ /* greyscale */ *r=v;*g=v;*b=v; return;}
  h /= 60; i=floor(h); /* six sectors */
  f = h-i;
  p = v*(1-s);
  q = v*(1-s*f);
  t = v*(1-s*(1-f));
  switch (i){
  case 0: *r=v;*g=t;*b=p; break;
  case 1: *r=q;*g=v;*b=p; break;
  case 2: *r=p;*g=v;*b=t; break;
  case 3: *r=p;*g=q;*b=v; break;
  case 4: *r=t;*g=p;*b=v; break;
  default: /* case 5 */ *r=v;*g=p;*b=q; break;}
}

void colorscale(double val,double valmin,double valmax,double *rcolor,double *gcolor,double *bcolor)
{
  /* use hsv with h in [240,-60] */
  double v=0;
  double h=0;
  if (valmax<=valmin){ valmax=valmin+1;}
  v = crop(val,valmin,valmax);
  h = 300*(valmax-v)/(valmax-valmin)-60;
  hsv2rgb(h,1,1,rcolor,gcolor,bcolor);
}

int WritePNMfile_color(double *ra,int rows,int cols,double min,double max,char *filename)
{
  /* This writes a double* array into a color pnm file at char* filename 
     assumes filename starts with "./" */
  int nr=0,nc=0;
  double d=0,rcolor=0,gcolor=0,bcolor=0;
  unsigned char *ud=NULL;
  FILE *fp;
  unsigned char *wkspace_mark=NULL;
  GLOBAL_wkspace_point = 0; wkspace_mark = wkspace_base;
  if (max<=min){ max=min+1;}
  if ((fp=fopen(filename,"w"))==NULL) { printf(" %% Warning: cannot open %s in WritePNMFile_color, writing to stdout\n",filename); fp=stdout;}
  fprintf(fp,"P6\n# %0.16lf %0.16lf\n%d %d\n%d\n",min,max,cols,rows,UCHAR_MAX);
  ud = (unsigned char *) wkspace_all0c(3*rows*cols);
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
      d = ra[nr+nc*rows]; colorscale(d,min,max,&rcolor,&gcolor,&bcolor);
      ud[0+nc*3+nr*3*cols] = (unsigned char)(UCHAR_MAX*rcolor); 
      ud[1+nc*3+nr*3*cols] = (unsigned char)(UCHAR_MAX*gcolor); 
      ud[2+nc*3+nr*3*cols] = (unsigned char)(UCHAR_MAX*bcolor);
      /* for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ }} */}}
  fwrite(ud,sizeof(unsigned char),3*rows*cols,fp);
  wkspace_reset(wkspace_mark); GLOBAL_wkspace_point = 1;
  if (fp!=stdout){ 
    fclose(fp);
    printf(" %% wrote file %s, (%d x %d pixels, %d bytes)\n",filename,rows,cols,(int)(rows*cols*sizeof(unsigned char)));
    return 1;}
  else{ printf(" %% wrote file %s to stdout\n",filename); return 0;}
  return 0;
}
