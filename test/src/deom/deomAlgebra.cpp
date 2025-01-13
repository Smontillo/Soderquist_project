#include "deom.hpp"

void deom::oprAct (cx_cube& d_ddos, const mat& sdip, const cube& pdip, const cx_vec& bdip, const cx_cube& ddos, const char lrc, const int xp) {

    const int nsav = nddo;
    d_ddos.slices(0,nddo-1).zeros();

	for (int iado=0; iado<nsav; ++iado) {

		const cx_mat& ado = ddos.slice(iado);

		if (iado==0 || is_valid (ado)) {

            hnod& nod = keys(iado);
            ivec key0 = nod.key;
			int loc0 = tree.find(key0)->rank;
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
                    ivec key1 = gen_key(key0, mp, 1);                 
                    const int m = modLabel(mp);
                    const int n = key1(mp)-1;
                    cx_double cl = sqrt((n+1)/coef_abs(mp))*coef_lft(mp);			
					//if(xp == 1) {d_ddos.slice(loc0) +=-deom_ci*(coef_lft(mp)-coef_rht(mp))*ado*qmd1.slice(m); }
                    if (!tree.try_insert(key1,nddo)) {	
                        int loc1 = tree.find(key1)->rank;
                        if (lrc == 'l') {
                            if (xp == 0) {d_ddos.slice(loc1) +=cl*pdip.slice(m)*ado;}
							else if (xp == 1) {d_ddos.slice(loc1) += cl*pdip.slice(m)*expn_gam(mp)*ado;}
                        } 
                    } else {
                        keys(nddo) = hnod(nod.gams+expn(mp),key1);
                        if (lrc == 'l'){
                            if (xp == 0) {d_ddos.slice(nddo) = cl*pdip.slice(m)*ado;}
							else if (xp == 1) {d_ddos.slice(nddo) =cl*pdip.slice(m)*expn_gam(mp)*ado;}
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
                        cx_double sn = sqrt(n*coef_abs(mp));
                        if (!tree.try_insert(key1,nddo)) {
                            int loc = tree.find(key1)->rank;
                            if (lrc == 'l'){
                                if (xp == 0) {d_ddos.slice(loc) += sn*pdip.slice(m)*ado;}
								else if (xp == 1){d_ddos.slice(loc) += -sn*expn_gam(mp)*pdip.slice(m)*ado;}
                            } 
                        } else {
                            keys(nddo) = hnod(nod.gams-expn(mp),key1);
                            if (lrc == 'l'){
                                if (xp == 0) {d_ddos.slice(nddo) = sn*pdip.slice(m)*ado;}
								else if (xp == 1) {d_ddos.slice(nddo) = -sn*expn_gam(mp)*pdip.slice(m)*ado;}
                            } 
                            nddo += 1;
                        }
                    }
                }
            }
        }
	}
}
