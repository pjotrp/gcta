/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for GWAS simulation
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

int gcta::read_QTL_file(string qtl_file, vector<string> &qtl_name, vector<int> &qtl_pos, vector<double> &qtl_eff, vector<int> &have_eff)
{
    qtl_name.clear();
    qtl_pos.clear();
    qtl_eff.clear();
    have_eff.clear();

    ifstream i_qtl(qtl_file.c_str());
    if(!i_qtl) throw("Error: can not open QTL file ["+qtl_file+"] to read.");
    string qtl_buf, str_buf;
    double qtl_eff_buf=0.0;
	cout<<"Reading a list of SNPs (as causal variants) from ["+qtl_file+"]."<<endl;
    map<string, int>::iterator iter, End=_snp_name_map.end();
    vector<string> vs_buf;
    vector<int> confirm(_snp_num);
    int icount=0;
	while(i_qtl){
		i_qtl>>qtl_buf;
		if(i_qtl.eof()) break;
		iter=_snp_name_map.find(qtl_buf);
		if( getline(i_qtl, str_buf) && StrFunc::split_string(str_buf, vs_buf)>0){
		    have_eff.push_back(1);
            qtl_eff_buf=atof(vs_buf[0].c_str());
            if(!(fabs(qtl_eff_buf)<1e5)) throw("Error: invalid effect size specified for the causal variant ["+str_buf+"].");
		}
		else{
		    have_eff.push_back(0);
		    qtl_eff_buf=0.0;
		}
		if(iter!=End){
		    qtl_name.push_back(qtl_buf);
		    qtl_pos.push_back(iter->second);
		    qtl_eff.push_back(qtl_eff_buf);
		}
	}
    i_qtl.close();
	cout<<qtl_pos.size()<<" SNPs (as causal variants) to be included from ["+qtl_file+"]."<<endl;
	return(qtl_pos.size());
}

void gcta::output_simu_par(vector<string> &qtl_name, vector<int> &qtl_pos, vector<double> &qtl_eff, double Vp)
{
    int i=0;
	vector<double> qsq(qtl_eff.size()), p(qtl_eff.size());
	for(i=0; i<qtl_eff.size(); i++){
        p[i]=0.5*_mu[qtl_pos[i]];
        qsq[i]=2.0*p[i]*(1.0-p[i])*qtl_eff[i]*qtl_eff[i];
        qsq[i]/=Vp;
    }
	string out_parfile=_out+".par";
	ofstream out_par(out_parfile.c_str());
	if(!out_par) throw("Error: can not open par file ["+out_parfile+"] to write!");
	out_par<<"QTL\tRefAllele\tFrequency\tEffect\tQsq"<<endl;
	for(i=0; i<qtl_eff.size(); i++) out_par<<qtl_name[i]<<"\t"<<_ref_A[qtl_pos[i]]<<"\t"<<p[i]<<"\t"<<qtl_eff[i]<<"\t"<<qsq[i]<<endl;
    out_par.close();
	cout<<"Simulated QTL effect(s) have been saved in ["+out_parfile+"]."<<endl;
}

void gcta::save_phenfile(vector< vector<double> > &y)
{
    string phenfile=_out+".phen";
	ofstream phen(phenfile.c_str());
	if(!phen) throw("Error: can not open the file ["+phenfile+"] to write.");
	int i=0, j=0;
	for(i=0; i<_keep.size(); i++){
		phen<<_fid[_keep[i]]<<" "<<_pid[_keep[i]]<<" ";
		for(j=0; j<y.size(); j++) phen<<y[j][i]<<" ";
		phen<<endl;
	}
	phen.close();
}

void gcta::GWAS_simu(string bfile, int simu_num, string qtl_file, int case_num, int control_num, double hsq, double K, bool output_causal, bool simu_emb_flag)
{
	int i=0, j=0;
	bool cc_flag=false;
	if(case_num>0 || control_num>0) cc_flag=true;

	cout<<"Simulation parameters:"<<endl;
	cout<<"Number of simulation replicate(s) = "<<simu_num<<" (Default = 1)"<<endl;
	cout<<"Heritability "<<(cc_flag?"of liability = ":" = ")<<hsq<<" (Default = 0.1)"<<endl;
	if(cc_flag){
	    cout<<"Disease prevalence = "<<K<<" (Default = 0.1)"<<endl;
	    cout<<"Number of cases = "<<case_num<<endl;
	    cout<<"Number of controls = "<<control_num<<endl;
	}
	cout<<endl;

    // Read QTL file
    vector<string> qtl_name;
    vector<int> qtl_pos, have_eff;
    vector<double> qtl_eff;
	int qtl_num=read_QTL_file(qtl_file, qtl_name, qtl_pos, qtl_eff, have_eff);
	update_id_map_kp(qtl_name, _snp_name_map, _include);

	// Generate QTL effects
	int Seed=-CommFunc::rand_seed();
	if(hsq>0.0){
		int num_gener_qtl_eff=0;
		for(i=0; i<qtl_num; i++){
			if(have_eff[i]==0){
				qtl_eff[i]=StatFunc::gasdev(Seed);
				num_gener_qtl_eff++;
			}
		}
		if(qtl_num-num_gener_qtl_eff>0) cout<<qtl_num-num_gener_qtl_eff<<" user-specified QTL effects."<<endl;
		if(num_gener_qtl_eff>0) cout<<num_gener_qtl_eff<<" unspecified QTL effects are generated from standard normal distribution."<<endl;
	}
	else{
		qtl_eff.clear();
		qtl_eff.resize(qtl_num);
	}


    // Calculate allele frequency
    vector< vector<float> > X;
    make_XMat(X, true);
    vector<double> sd_SNP_i;
    std_XMat(X, sd_SNP_i, false, 0);

	// Calculate Ve and threhold
    double var_g=0.0, var_e=1.0;
    vector<double> g(_keep.size());
    if(hsq>0.0){
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<qtl_num; j++) g[i]+=X[i][j]*qtl_eff[j];
        }
        var_g=CommFunc::var(g);
	    var_e=var_g*(1.0/hsq-1.0);
    }
	double sd_e=sqrt(var_e);

    // Output par file
    output_simu_par(qtl_name, qtl_pos, qtl_eff, var_e+var_g);

    // Output phenotype file
    cout<<"Simulating GWAS based on the real genotyped data with "<<simu_num<<" replicate(s) ..."<<endl;
    vector< vector<double> > y(simu_num);
    int case_num_buf=0, control_num_buf=0;
    for(i=0; i<simu_num; i++){
        y[i].resize(_keep.size());
        for(j=0; j<_keep.size(); j++) y[i][j]=g[j]+sd_e*StatFunc::gasdev(Seed);
        if(cc_flag){
            case_num_buf=0; control_num_buf=0;
            vector<double> y_buf(y[i]);
            stable_sort(y_buf.begin(), y_buf.end());
            int n=(int)(_indi_num*(1.0-K));
            double Th=0.5*(y_buf[n]+y_buf[n-1]);
            for(j=0; j<_keep.size(); j++){
                if(y[i][j]>Th){
                    if(case_num_buf<case_num){ y[i][j]=2; case_num_buf++; }
                    else y[i][j]=-9;
                }
                else{
                    if(control_num_buf<control_num){ y[i][j]=1; control_num_buf++; }
                    else y[i][j]=-9;
                }
            }
        }
    }

	if(!simu_emb_flag){
	    save_phenfile(y);
	    if(cc_flag) cout<<"Simulated "<<case_num_buf<<" cases and "<<control_num<<" controls have been saved in ["+_out+".phen"+"]."<<endl;
	    else cout<<"Simulated phenotypes of "<<_keep.size()<<" individuals have been saved in ["+_out+".phen"+"]."<<endl;
	}
	else{
        // emBayesB format
        if(!output_causal) update_id_map_rm(qtl_name, _snp_name_map, _include);
        string out_rstfile=_out+".emb";
        ofstream out_emBayesB(out_rstfile.c_str());
        if(!out_emBayesB) throw("Error: can not open the file ["+out_rstfile+"] to write.");
        cout<<"Saving the simulated data to the file ["+out_rstfile+"] (in emBayesB format)."<<endl;
        for(i=0; i<_keep.size(); i++){
            if(y[0][i]==-9) continue;
            out_emBayesB<<_pid[_keep[i]]<<" "<<g[i]<<" "<<y[0][i]<<endl;
            for(j=0; j<_include.size(); j++){
                if(_snp_1[_include[j]][i] && !_snp_2[_include[j]][i]) out_emBayesB<<_mu[_include[j]]<<" ";
                else out_emBayesB<<(double)(_snp_1[_include[j]][i]+_snp_2[_include[j]][i])<<" ";
            }
            out_emBayesB<<endl;
        }
        out_emBayesB.close();
        cout<<"Simulated data ("<<_keep.size()<<" individuals and "<<_include.size()<<" SNPs) has been saved in ["+out_rstfile+"]."<<endl;
	}
}

void gcta::genet_dst(string bfile, string hapmap_genet_map)
{
	// Read bim file
	read_bimfile(bfile+".bim");
    int snp_num=_snp_name.size();

    // Read HAPMAP genetic map files
    int i=0, j=0;
	string str_buf;
	vector<string> vs_buf;
    string genet_mapfile;
    vector< vector< vector<double> > > hap_genet(_autosome_num);
    for(i=0; i<_autosome_num; i++){
        stringstream str_strm;
        str_strm<<hapmap_genet_map<<i+1<<"_CEU_b36.txt";
        genet_mapfile=str_strm.str();
        ifstream i_genet_map(genet_mapfile.c_str());
        if(!i_genet_map) throw("Error: can not open HAPMAP genetic map file "+genet_mapfile+"!");
        hap_genet[i].resize(2);
        getline(i_genet_map, str_buf);
    	while(getline(i_genet_map, str_buf)){
    		if(StrFunc::split_string(str_buf, vs_buf) < 3) continue;
    		hap_genet[i][0].push_back(atof(vs_buf[0].c_str()));
    		hap_genet[i][1].push_back(atof(vs_buf[1].c_str()));
    	}
    	i_genet_map.close();
    }

    // calculate genetic distance
    int pos1=0, pos2=0;
    double prev_bp=0.0;
    vector<double> dst(snp_num);
    vector<double>::iterator iter1, iter2;
    for(i=0; i<snp_num; i++){
       if(i==0 || _chr[i-1]!=_chr[i]){ dst[i]=0.0; continue; }
       iter1=upper_bound(hap_genet[_chr[i]-1][0].begin(), hap_genet[_chr[i]-1][0].end(), _bp[i-1]);
       iter2=upper_bound(hap_genet[_chr[i]-1][0].begin(), hap_genet[_chr[i]-1][0].end(), _bp[i]);
       if(iter1==hap_genet[_chr[i]-1][0].end() || iter2==hap_genet[_chr[i]-1][0].end()) { dst[i]=_bp[i]-_bp[i-1]; continue; }
       pos1=iter1-hap_genet[_chr[i]-1][0].begin();
       pos2=iter2-hap_genet[_chr[i]-1][0].begin();
       prev_bp=_bp[i-1];
       while(pos1<pos2){
           dst[i]+=(hap_genet[_chr[i]-1][0][pos1]-prev_bp)*hap_genet[_chr[i]-1][1][pos1-1];
           prev_bp=hap_genet[_chr[i]-1][0][pos1];
           pos1++;
       }
       dst[i]+=(_bp[i]-prev_bp)*hap_genet[_chr[i]-1][1][pos1-1];
    }

    // Output fam file
	string out_bimfile=_out+".genetdst";
	ofstream out_bim(out_bimfile.c_str());
	if(!out_bim) throw("Error: can not open file "+out_bimfile+" to write!");
	for(i=0; i<snp_num; i++) out_bim<<_chr[i]<<"\t"<<_snp_name[i]<<"\t"<<dst[i]*1e-6<<"\t"<<_bp[i]<<"\t"<<_allele1[i]<<"\t"<<_allele2[i]<<endl;
	out_bim.close();
	cout<<"Genetic distances have been created, and been saved in ["+out_bimfile+"]."<<endl;
}
