#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

unsigned char *wkspace_alloc_nocheck(unsigned long long int size) {
  unsigned char* retval;
  if (wkspace_left < size) { printf(" %% Warning! out of memory (%lld < %llu)\n",wkspace_left,size); exit(RET_NOMEM); return NULL;}
  size = CACHEALIGN(size); retval = wkspace_base; wkspace_base += size; wkspace_left -= size; wkspace_used += size;
  return retval;
}

struct wkspace_point * wkspace_make_point()
{
  struct wkspace_point *w=NULL;
  w = (struct wkspace_point *)wkspace_alloc_nocheck(sizeof(struct wkspace_point)*1);
  w->parent = NULL;
  w->child = NULL;
  w->point = (long long int *)wkspace_alloc_nocheck(sizeof(long long int)*1);
  w->check = 0; *(w->point) = w->check;
  return w;
}

struct wkspace_point * wkspace_set_point(struct wkspace_point *w)
{
  w->child = (struct wkspace_point *)wkspace_alloc_nocheck(sizeof(struct wkspace_point)*1);
  w->child->parent = w;
  w->child->child = NULL;
  w->child->point = (long long int *)wkspace_alloc_nocheck(sizeof(long long int)*1);
  w->child->check = w->check + 1; *(w->child->point) = w->child->check;
  return w->child;
}

void wkspace_printf_point(struct wkspace_point *w){
  if (GLOBAL_wkspace_point){
    printf(" %% parent %p, self %p, child %p, check %lld point %lld\n",w->parent,w,w->child,w->check,*(w->point));
    if (w->child!=NULL){ wkspace_printf_point(w->child);}
    /* if (GLOBAL_wkspace_point){ } */}
}

long long int wkspace_check_point(struct wkspace_point *w){
  if (GLOBAL_wkspace_point){  
    if (w->check == *(w->point)){ if (w->child!=NULL){ return wkspace_check_point(w->child);} else{ return w->check;}}
    else{ printf(" %% Warning! parent %p, self %p, child %p, check %lld, point %lld\n",w->parent,w,w->child,w->check,*(w->point)); exit(1); return 0;}
    /* if (GLOBAL_wkspace_point){ } */}
  return 0;
}

unsigned char* wkspace_alloc(unsigned long long int size) {
  int verbose=0;
  unsigned char* retval;
  if (wkspace_left < size) { printf(" %% Warning! out of memory (%lld < %llu)\n",wkspace_left,size); exit(RET_NOMEM); return NULL;}
  if (verbose){ printf(" %% [entering wkspace_alloc] wkspace_alloc: wkspace_base %p; mod 16: %s\n",wkspace_base,((uintptr_t)wkspace_base/16)*16==(uintptr_t)wkspace_base ? "Yes" : "No!");}
  if (verbose){ printf(" %% wkspace_alloc: size %llu --> ",size);}
  size = CACHEALIGN(size); if (verbose){ printf(" size %llu\n",size);}
  retval = wkspace_base; wkspace_base += size; wkspace_left -= size; wkspace_used += size;
  if (verbose){ printf(" %% [finished wkspace_alloc] wkspace_alloc: wkspace_base %p; mod 16: %s\n",wkspace_base,((uintptr_t)wkspace_base/16)*16==(uintptr_t)wkspace_base ? "Yes" : "No!");}
  if (GLOBAL_wkspace_point){ wkspace_point_t = wkspace_set_point(wkspace_point_t);}
  return retval;
}

unsigned char* wkspace_all0c(unsigned long long int size) {
  int verbose=0;
  unsigned char* retval;
  if (wkspace_left < size) { printf(" %% Warning! out of memory (%lld < %llu)\n",wkspace_left,size); exit(RET_NOMEM); return NULL;}
  if (verbose){ printf(" %% [entering wkspace_alloc] wkspace_alloc: wkspace_base %p; mod 16: %s\n",wkspace_base,((uintptr_t)wkspace_base/16)*16==(uintptr_t)wkspace_base ? "Yes" : "No!");}
  if (verbose){ printf(" %% wkspace_alloc: size %llu --> ",size);}
  size = CACHEALIGN(size); if (verbose){ printf(" size %llu\n",size);}
  retval = wkspace_base; wkspace_base += size; wkspace_left -= size; wkspace_used += size;
  if (verbose){ printf(" %% [finished wkspace_alloc] wkspace_alloc: wkspace_base %p; mod 16: %s\n",wkspace_base,((uintptr_t)wkspace_base/16)*16==(uintptr_t)wkspace_base ? "Yes" : "No!");}
  fill_uchar_zero(retval,size);
  if (GLOBAL_wkspace_point){ wkspace_point_t = wkspace_set_point(wkspace_point_t);}
  return retval;
}

void wkspace_reset(void* new_base) {
  unsigned long freed_bytes = wkspace_base - (unsigned char*)new_base;
  wkspace_base = (unsigned char*)new_base;
  wkspace_left += freed_bytes; wkspace_used -= freed_bytes;
}

void wkspace_printf(){ 
  printf(" %% CHECK:%lld used/left (%lld/%lld UB) (%.2f/%.2f KB) (%.2f/%.2f MB) (%.2f/%.2f GB)\n",wkspace_check_point(wkspace_point_0),wkspace_used,wkspace_left,(double)wkspace_used/1024,(double)wkspace_left/1024,(double)wkspace_used/1048576,(double)wkspace_left/1048576,(double)wkspace_used/1073741824,(double)wkspace_left/1073741824);
}
