#include "deom.hpp"

void deom::oprAct (cx_cube& d_ddos, const mat& sdip, const cube& pdip, const vec& bdip, const cx_cube& ddos, const char lrc, const int xp) {

    const int nsav = nddo;
    d_ddos.slices(0,nddo-1).zeros();

	const int nper=nind/nmod;

    for (int iado=0; iado<nsav; ++iado) {

        const cx_mat& ado = ddos.slice(iado);

        if (iado==0 || is_valid (ado)) {

            hnod& nod = keys(iado);
            ivec key0 = nod.key;
            int tier = tree.find(key0)->tier;

            if (lrc == 'l') {
                d_ddos.slice(iado) += sdip*ado;
            } else if (lrc == 'r') {
                d_ddos.slice(iado) += ado*sdip;
            } else if (lrc == 'c') {
                d_ddos.slice(iado) += sdip*ado-ado*sdip;
            }

            if (tier < lmax) {
                for (int mp=0; mp<nind; ++mp) {
                    int mpp=0;
					if ((mp%nper)==0){mpp=mp+1;}
					else if ((mp%nper)==1) {mpp=mp-1;}
					else {mpp=mp;}
                    ivec key1 = gen_key(key0, mpp, 1);
                    const int m = modLabel(mp);
                    const int n = key1(mpp)-1;
                    cx_double cl = sqrt((n+1)/coef_abs(mpp))*coef_lft(mpp)*bdip(mp);			
                    cx_double cr = sqrt((n+1)/coef_abs(mpp))*coef_rht(mpp)*bdip(mp);
					cx_double cll = sqrt((n+1)/coef_abs(mpp))*coef_lft(mp)*bdip(mp);				
					cx_double crr = sqrt((n+1)/coef_abs(mpp))*coef_rht(mp)*bdip(mp);
                    if (!tree.try_insert(key1,nddo)) {
                        int loc = tree.find(key1)->rank;
                        if (lrc == 'l') {
                            if (xp == 0) {d_ddos.slice(loc) += cl*pdip.slice(m)*ado;}
							else if (xp == 1) {d_ddos.slice(loc) += -crr*expn_gam(mp)*pdip.slice(m)*ado;}
							//else if (xp == 1&& ((mp%nper)!=0 && (mp%nper)!=1)) {d_ddos.slice(loc) += deom_ci*sqrt((n+1)/coef_abs(mp))*ado*coef_lft(mp);}
                        } else if (lrc == 'r') {
                            if (xp == 0) {d_ddos.slice(loc) += cr*ado*pdip.slice(m);}
							//else if ((xp == 1) && ((mp%nper)==0 || (mp%nper)==1)) {d_ddos.slice(loc) += -cll*expn_gam(mp)*ado*pdip.slice(m);}
                        } else if (lrc == 'c') {
                            if (xp == 0) {d_ddos.slice(loc) += cl*pdip.slice(m)*ado-cr*ado*pdip.slice(m);}
							//else if ((xp == 1) && ((mp%nper)==0 || (mp%nper)==1)){d_ddos.slice(loc) += -crr*expn_gam(mp)*pdip.slice(m)*ado+cll*expn_gam(mp)*ado*pdip.slice(m);}
                        }
                    } else {
                        keys(nddo) = hnod(nod.gams+expn(mpp),key1);
                        if (lrc == 'l'){
                            if (xp == 0) {d_ddos.slice(nddo) = cl*pdip.slice(m)*ado;}
							else if (xp == 1) {d_ddos.slice(nddo) =-crr*expn_gam(mp)*pdip.slice(m)*ado;}
							//else if (xp == 1&& ((mp%nper)!=0 && (mp%nper)!=1)) {d_ddos.slice(nddo) = deom_ci*sqrt((n+1)/coef_abs(mp))*ado*coef_lft(mp);}
                        } else if (lrc == 'r') {
                            if (xp == 0) {d_ddos.slice(nddo) = cr*ado*pdip.slice(m);}
							//else if ((xp == 1) && ((mp%nper)==0 || (mp%nper)==1)) {d_ddos.slice(nddo) = -cll*expn_gam(mp)*ado*pdip.slice(m);}
                        } else if (lrc == 'c') {
                            if (xp == 0) {d_ddos.slice(nddo) = cl*pdip.slice(m)*ado-cr*ado*pdip.slice(m);}
							//else if ((xp == 1) && ((mp%nper)==0 || (mp%nper)==1)) {d_ddos.slice(nddo) = -crr*expn_gam(mp)*pdip.slice(m)*ado+cll*expn_gam(mp)*ado*pdip.slice(m);}
                        }
                        nddo += 1;
                    }
                }
            }

            if (iado > 0) {
                for (int mp=0; mp<nind; ++mp) {
                    ivec key1 = gen_key(key0, mp, -1);
                    if (!key1.is_empty()) {
                        const int m = modLabel(mp);
                        const int n = (key1.n_rows<(unsigned int)(mp+1))?1:(key1(mp)+1);
                        cx_double sn = sqrt(n*coef_abs(mp))*bdip(mp);
                        if (!tree.try_insert(key1,nddo)) {
                            int loc = tree.find(key1)->rank;
                            if (lrc == 'l'){
                                if (xp == 0) {d_ddos.slice(loc) += sn*pdip.slice(m)*ado;}
								else if (xp == 1){d_ddos.slice(loc) += -sn*expn_gam(mp)*pdip.slice(m)*ado;}
                            } else if (lrc == 'r') {
                                if (xp == 0) {d_ddos.slice(loc) += sn*ado*pdip.slice(m);}
								//else if ((xp == 1) && ((mp%nper)==0 || (mp%nper)==1)){d_ddos.slice(loc) += -sn*expn_gam(mp)*ado*pdip.slice(m);}
                            } else if (lrc == 'c'){ 
                                if (xp == 0) {d_ddos.slice(loc) += sn*(pdip.slice(m)*ado-ado*pdip.slice(m));}
								//else if ((xp == 1) && ((mp%nper)==0 || (mp%nper)==1)) {d_ddos.slice(loc) += -sn*expn_gam(mp)*(pdip.slice(m)*ado-ado*pdip.slice(m));}
                            }
                        } else {
                            keys(nddo) = hnod(nod.gams-expn(mp),key1);
                            if (lrc == 'l'){
                                if (xp == 0) {d_ddos.slice(nddo) = sn*pdip.slice(m)*ado;}
								else if (xp == 1) {d_ddos.slice(nddo) = -sn*expn_gam(mp)*pdip.slice(m)*ado;}
                            } else if (lrc == 'r') {
                                if (xp == 0) {d_ddos.slice(nddo) = sn*ado*pdip.slice(m);}
								//else if ((xp == 1) && ((mp%nper)==0 || (mp%nper)==1)) {d_ddos.slice(nddo) = -sn*expn_gam(mp)*ado*pdip.slice(m);}
                            } else if (lrc == 'c') {
                                if (xp == 0) {d_ddos.slice(nddo) = sn*(pdip.slice(m)*ado-ado*pdip.slice(m));}
								//else if ((xp == 1) && ((mp%nper)==0 || (mp%nper)==1)) {d_ddos.slice(nddo) = -sn*expn_gam(mp)*(pdip.slice(m)*ado-ado*pdip.slice(m));}
                            }
                            nddo += 1;
                        }
                    }
                }
            }
        }
    }
}
