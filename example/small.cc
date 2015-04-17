#include "opm.h"
#include <algorithm>
#include <stdlib.h>
#include <signal.h>

#include "data/small_ref.h"
#include "data/medium_ref.h"
#include "data/large_ref.h"

#include "data/small_permute.h"
#include "data/small_posrotate.h"
#include "data/small_posshift.h"
#include "data/small_posrotshift.h"
#include "data/small_target.h"


FILE *fd = NULL;
wcsprm* wcs;
wcsprm* new_wcs;

void
killgnuplot(int s)
{
  FILE* pid = popen("ps | awk '/gnuplot/{print $1}'","r");
  int id;
  if (fscanf(pid,"%d",&id)>0)
    kill(id, SIGTERM);
  exit(1);
}

void
run_matching
(const double Nmatched,
 const std::vector<opm::xym> &ref_,
 const std::vector<opm::xym> &obj_,
 const char* filename)
{
  wcsini(0, 2, wcs);
  wcs->crpix[0] = 512.0;
  wcs->crpix[1] = 512.0;
  wcs->cdelt[0] = -0.00026278;
  wcs->cdelt[1] =  0.00026278;
  wcs->crval[0] =  160.;
  wcs->crval[1] =   35.;
  wcs->equinox  = 2000.;
  sprintf(wcs->cunit[0],"deg");
  sprintf(wcs->cunit[1],"deg");
  sprintf(wcs->ctype[0],"RA---TAN");
  sprintf(wcs->ctype[1],"DEC--TAN");
  sprintf(wcs->wcsname,"FK5");
  wcsset(wcs);

  opm::referencelist ref
    = opm::generate_referencelist(ref_,wcs);
  opm::database refdb = opm::generate_database(ref);

  opm::objectlist obj
    = opm::generate_objectlist(obj_);
  opm::triangles objdb = opm::generate_triangles(obj);

  opm::matched_list matched;
  for (auto &p : objdb) {
    std::vector<opm::conversion> x = opm::pre::match(p, .1, .1, refdb);
    for (auto &v : x) {
      matched = opm::final::match(v, obj, 5., ref);
      fprintf(stderr,"matched.size(): %ld\n", matched.size());
      if (matched.size() >= Nmatched) goto match_success;
    }
  }
  delete(new_wcs);
  delete(wcs);
  pclose(fd);
  killgnuplot(1);
  exit(1);

 match_success:
  printf("## CROSS CHECK ##\n");
  wcscopy(1, wcs, new_wcs);
  update_wcsprm(new_wcs, matched);
  fprintf(fd,"$DATA << EOD\n");
  for (auto &p : matched) {
    int status;
    double world[2], new_world[2];
    double imgcrd[2], phi[1], theta[1];
    double pixcrd[2] = {p.first.x, p.first.y};
    wcsp2s(wcs,     1, 2, pixcrd, imgcrd, phi, theta, world, &status);
    wcsp2s(new_wcs, 1, 2, pixcrd, imgcrd, phi, theta, new_world, &status);
    printf("(%12.8lf %12.8lf) (%12.8lf %12.8lf) [%+lf %+lf]\n",
           world[0],world[1],new_world[0],new_world[1],
           3600.*(new_world[0]-p.second.x),3600.*(new_world[1]-p.second.y));
    fprintf(fd,"%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n",
            world[0],world[1], new_world[0], new_world[1],
            p.second.x, p.second.y);
  }
  printf("## CROSS CHECK ##\n");
  fprintf(fd,"EOD\n" );
  fprintf(fd, "plot '%s' u 1:2 not ","data/small_ref.sav");
  fprintf(fd, "w p lw 0.2 ps 2 pt 6\n");
  fprintf(fd, "replot '%s' u 1:2 t 'REFERENCE' ", filename);
  fprintf(fd, "w p ps 0.5 pt 5; pause 1\n");
  fprintf(fd, "replot $DATA u 1:2 t 'OBJECT' ");
  fprintf(fd, "w p ps 2 lw 2 pt 1; pause 1\n");
  fprintf(fd, "replot $DATA u 1:2:($3-$1):($4-$2) t 'Match' w vector ");
  fprintf(fd, "filled lw 2 \n");
  fprintf(fd,"pause 3\n");
  fflush(fd);
}

int
main(int argc, char *argv[])
{
  wcs = new wcsprm;
  new_wcs = new wcsprm;
  fd = popen("gnuplot","w");
  signal(SIGINT, killgnuplot);

  fprintf(fd,"set terminal wxt enh font ',12'\n");
  fprintf(fd,"set key rmargin vertical Left\n");
  fprintf(fd,"set key samplen 1\n");
  fprintf(fd,"set grid x y mx my\n");
  fprintf(fd,"set grid lt 0 lw 1, lt 0 lw 0.5\n");
  fprintf(fd,"set xr [157:163]\n");
  fprintf(fd,"set yr [33:37]\n");
  fprintf(fd,"set xtics 0,1,200 format '%%.0f'\n");
  fprintf(fd,"set ytics 0,1,200 format '%%.0f'\n");
  fprintf(fd,"set mxtics 5\n");
  fprintf(fd,"set mytics 5\n");
  
  run_matching(6, opm_test::small_ref, opm_test::small_permute,
               "data/small_ref.sav"); exit(1);
  run_matching(6, opm_test::small_ref, opm_test::small_posshift,
               "data/small_ref.sav");
  run_matching(6, opm_test::small_ref, opm_test::small_posrotate,
               "data/small_ref.sav");
  run_matching(6, opm_test::small_ref, opm_test::small_posrotshift,
               "data/small_ref.sav");
  run_matching(6, opm_test::small_ref, opm_test::small_target,
               "data/small_ref.sav");

  run_matching(6, opm_test::medium_ref, opm_test::small_permute,
               "data/medium_ref.sav");
  run_matching(6, opm_test::medium_ref, opm_test::small_posshift,
               "data/medium_ref.sav");
  run_matching(6, opm_test::medium_ref, opm_test::small_posrotate,
               "data/medium_ref.sav");
  run_matching(6, opm_test::medium_ref, opm_test::small_posrotshift,
               "data/medium_ref.sav");
  run_matching(6, opm_test::medium_ref, opm_test::small_target,
               "data/medium_ref.sav");

  run_matching(6, opm_test::large_ref, opm_test::small_permute,
               "data/large_ref.sav");
  run_matching(6, opm_test::large_ref, opm_test::small_posshift,
               "data/large_ref.sav");
  run_matching(6, opm_test::large_ref, opm_test::small_posrotate,
               "data/large_ref.sav");
  run_matching(6, opm_test::large_ref, opm_test::small_posrotshift,
               "data/large_ref.sav");
  run_matching(6, opm_test::large_ref, opm_test::small_target,
               "data/large_ref.sav");

  delete(new_wcs);
  delete(wcs);
  pclose(fd);
  return 0;
}
