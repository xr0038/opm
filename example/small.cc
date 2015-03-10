#include "opm.h"
#include <algorithm>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>

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

void
run_matching
(const double Nmatched,
 const std::vector<opm::xym> &ref_,
 const std::vector<opm::xym> &obj_,
 const char* filename_1, const char* filename_2)
{
  opm::referencelist ref
    = opm::generate_referencelist(ref_,wcs);
  opm::database refdb = opm::generate_database(ref);

  opm::objectlist obj
    = opm::generate_objectlist(obj_);
  opm::triangles objdb = opm::generate_triangles(obj);

  opm::matched_list matched;
  for (auto &p : objdb) {
    std::vector<opm::conversion> x = opm::pre::match(p, .3, .3, refdb);
    for (auto &v : x) {
      fprintf(stderr,"%.5lf,%.5lf,%.5lf,%.5lf,%.5lf,%.5lf\n",
              v[0],v[1],v[2],v[3],v[4],v[5]);
      matched = opm::final::match(v, obj, 1.0e2, ref);
      fprintf(stderr,"matched.size(): %ld\n", matched.size());
      if (matched.size() >= Nmatched) goto match_success;
    }
  }
  return;

 match_success:
  fprintf(fd, "plot '%s' u 1:2 not ","data/small_permute.sav");
  fprintf(fd, "w p lw 0.2 ps 2 pt 6\n");
  fprintf(fd, "replot '%s' u 1:2 t 'REFERENCE' ", filename_1);
  fprintf(fd, "w p ps 1 pt 5; pause 1\n");
  fprintf(fd, "replot '%s' u 1:2 t 'OBJECT' ", filename_2);
  fprintf(fd, "w p ps 2 lw 2 pt 1; pause 1\n");
  fprintf(fd, "replot '-' u 1:2:($3-$1):($4-$2) t 'Match' w vector ");
  fprintf(fd, "filled lw 2 \n");
  for (auto &p : matched) {
    fprintf(fd,"%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n",
            p.first.x, p.first.y,
            p.second.xi(), p.second.eta(),
            p.second.x, p.second.y);
  }
  fprintf(fd,"e\n" );
  fprintf(fd,"pause 3\n");
  fflush(fd);
}

void
killgnuplot(int s)
{
  FILE* pid = popen("ps | awk '/gnuplot/{print $1}'","r");
  int id;
  if (fscanf(pid,"%d",&id)>0)
    kill(id, SIGTERM);
  exit(1);
}

int
main(int argc, char *argv[])
{
  wcs = new wcsprm;
  wcsini(0, 2, wcs);
  wcs->crpix[0] = 512;
  wcs->crpix[1] = 512;
  wcs->cdelt[0] = -0.00026278;
  wcs->cdelt[1] =  0.00026278;
  wcs->crval[0] =  160.;
  wcs->crval[1] =   35.;
  sprintf(wcs->cunit[0],"deg");
  sprintf(wcs->cunit[1],"deg");
  sprintf(wcs->ctype[0],"RA---TAN");
  sprintf(wcs->ctype[1],"DEC--TAN");
  sprintf(wcs->wcsname,"FK5");
  wcsset(wcs);


  fd = popen("gnuplot","w");
  signal(SIGINT, killgnuplot);

  fprintf(fd,"set terminal wxt enh font ',12'\n");
  fprintf(fd,"set key rmargin vertical Left\n");
  fprintf(fd,"set key samplen 1\n");
  fprintf(fd,"set grid x y mx my\n");
  fprintf(fd,"set xr [-7000:7000]\n");
  fprintf(fd,"set yr [-7000:7000]\n");
  
  run_matching(6, opm_test::small_ref, opm_test::small_permute,
               "data/small_permute.sav", "data/small_permute.sav");
  run_matching(6, opm_test::small_ref, opm_test::small_posshift,
               "data/small_permute.sav", "data/small_posshift.sav");
  run_matching(6, opm_test::small_ref, opm_test::small_posrotate,
               "data/small_permute.sav", "data/small_posrotate.sav");
  run_matching(6, opm_test::small_ref, opm_test::small_posrotshift,
               "data/small_permute.sav", "data/small_posrotshift.sav");
  run_matching(6, opm_test::small_ref, opm_test::small_target,
               "data/small_permute.sav", "data/small_target.sav");

  run_matching(6, opm_test::medium_ref, opm_test::small_permute,
               "data/medium_permute.sav", "data/small_permute.sav");
  run_matching(6, opm_test::medium_ref, opm_test::small_posshift,
               "data/medium_permute.sav", "data/small_posshift.sav");
  run_matching(6, opm_test::medium_ref, opm_test::small_posrotate,
               "data/medium_permute.sav", "data/small_posrotate.sav");
  run_matching(6, opm_test::medium_ref, opm_test::small_posrotshift,
               "data/medium_permute.sav", "data/small_posrotshift.sav");
  run_matching(6, opm_test::medium_ref, opm_test::small_target,
               "data/medium_permute.sav", "data/small_target.sav");

  run_matching(6, opm_test::large_ref, opm_test::small_permute,
               "data/large_permute.sav", "data/small_permute.sav");
  run_matching(6, opm_test::large_ref, opm_test::small_posshift,
               "data/large_permute.sav", "data/small_posshift.sav");
  run_matching(6, opm_test::large_ref, opm_test::small_posrotate,
               "data/large_permute.sav", "data/small_posrotate.sav");
  run_matching(6, opm_test::large_ref, opm_test::small_posrotshift,
               "data/large_permute.sav", "data/small_posrotshift.sav");
  run_matching(6, opm_test::large_ref, opm_test::small_target,
               "data/large_permute.sav", "data/small_target.sav");

  pclose(fd);
  delete wcs;
  return 0;
}
