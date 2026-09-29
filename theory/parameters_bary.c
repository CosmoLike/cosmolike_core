void set_bary_parameters_to_mb2();
void set_bary_parameters_to_illustris();
void set_bary_parameters_to_HzAGN();
void set_bary_parameters_to_eagle();

void set_bary_parameters_to_owls_AGN();
void set_bary_parameters_to_owls_DBLIMFV1618();
void set_bary_parameters_to_owls_NOSN();
void set_bary_parameters_to_owls_NOSN_NOZCOOL();
void set_bary_parameters_to_owls_NOZCOOL();
void set_bary_parameters_to_owls_REF();
void set_bary_parameters_to_owls_WDENS();
void set_bary_parameters_to_owls_WML1V848();
void set_bary_parameters_to_owls_WML4();
void set_bary_parameters_to_antilles(int id);
void reset_PkRatio_baryons(void);

void set_bary_parameters_to_TNG100()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_TNG100.dat");
	sprintf(bary.scenario,"TNG100");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 13;		
	double z[13] = {3.71,3.49,3.28,2.90,2.44,2.1,1.74,1.41,1.04,0.7,0.35,0.18,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}


void set_bary_parameters_to_HzAGN()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_HzAGN.dat");
	sprintf(bary.scenario,"HzAGN");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 11;		
	double z[11] = {4.9285,4.249,3.7384,3.33445,3.00295,1.96615,1.02715,0.519195,0.22878,0.017865,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}


void set_bary_parameters_to_mb2()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_mb2.dat");
	sprintf(bary.scenario,"mb2");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 21;		
	double z[21] = {3.5,3.25,2.8,2.45,2.1,2.0,1.8,1.7,1.6,1.4,1.2,1.1,1.0,0.8,0.7,0.6,0.4,0.35,0.2,0.0625,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_illustris()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_ill1.dat");
	sprintf(bary.scenario,"illustris");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 23;		
	double z[23] = {3.5,3.49,3.28,3.08,2.90,2.73,2.44,2.1,2.0,1.82,1.74,1.6,1.41,1.21,1.04,1.0,0.79,0.7,0.6,0.4,0.35,0.2,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_eagle()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_eagle.dat");
	sprintf(bary.scenario,"eagle");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 13;		
	double z[13] = {3.53,3.02,2.48,2.24,2.01,1.74,1.49,1.26,1.0,0.74,0.5,0.27,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_owls_AGN()
{
	//sprintf(bary.FILE_logPkR,"./theory_hj/data_sims/PkRatio_0.3_10/PkRatio_owls_AGN.dat");
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_owls_AGN.dat");
	//sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_owls_AGN.dat");
	sprintf(bary.scenario,"owls_AGN");
	//sprintf(bary.FILE_logPkR,"./theory_hj/data_sims/PkRatio_owls_AGN.dat");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 16;		
	double z[16] = {3.5,3.25,3.0,2.75,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_owls_DBLIMFV1618()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_owls_DBLIMFV1618.dat");
	sprintf(bary.scenario,"owls_DBLIMFV1618");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 16;		
	double z[16] = {3.5,3.25,3.0,2.75,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_owls_NOSN()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_owls_NOSN.dat");
	sprintf(bary.scenario,"owls_NOSN");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 15;		
	double z[15] = {3.5,3.25,3.0,2.75,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.0};   //no z=0.125 col for NOSN
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_owls_NOSN_NOZCOOL()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_owls_NOSN_NOZCOOL.dat");
	sprintf(bary.scenario,"owls_NOSN_NOZCOOL");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 16;		
	double z[16] = {3.5,3.25,3.0,2.75,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_owls_NOZCOOL()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_owls_NOZCOOL.dat");
	sprintf(bary.scenario,"owls_NOZCOOL");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 16;		
	double z[16] = {3.5,3.25,3.0,2.75,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_owls_REF()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_owls_REF.dat");
	sprintf(bary.scenario,"owls_REF");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 16;		
	double z[16] = {3.5,3.25,3.0,2.75,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_owls_WDENS()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_owls_WDENS.dat");
	sprintf(bary.scenario,"owls_WDENS");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 16;		
	double z[16] = {3.5,3.25,3.0,2.75,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_owls_WML1V848()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_owls_WML1V848.dat");
	sprintf(bary.scenario,"owls_WML1V848");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 16;		
	double z[16] = {3.5,3.25,3.0,2.75,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_owls_WML4()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_owls_WML4.dat");
	sprintf(bary.scenario,"owls_WML4");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 16;		
	double z[16] = {3.5,3.25,3.0,2.75,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
	
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_cowls_AGN()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_cowls_AGN_Mseed800_WMAP7_L100N512.dat");
	sprintf(bary.scenario,"cowls_AGN");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 15;		
	double z[15] = {3.0,2.75,2.50,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
		
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_cowls_AGN_T8p5()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_cowls_AGN_Mseed800_T8.5_WMAP7_L100N512.dat");
	sprintf(bary.scenario,"cowls_AGN_T8p5");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 15;		
	double z[15] = {3.0,2.75,2.50,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
		
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_cowls_AGN_T8p7()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_cowls_AGN_Mseed800_T8.7_WMAP7_L100N512.dat");
	sprintf(bary.scenario,"cowls_AGN_T8p7");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 15;		
	double z[15] = {3.0,2.75,2.50,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
		
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_BAHAMAS()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_BAHAMAS_nu0_WMAP9_L400N1024.dat");
	sprintf(bary.scenario,"BAHAMAS");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 15;		
	double z[15] = {3.0,2.75,2.50,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
		
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_BAHAMAS_T7p6()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_BAHAMAS_T7.6_WMAP9_L400N1024.dat");
	sprintf(bary.scenario,"BAHAMAS_T7p6");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 15;		
	double z[15] = {3.0,2.75,2.50,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
		
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_BAHAMAS_T8p0()
{
	sprintf(bary.FILE_logPkR,"../cosmolike_core/logPkRatio/logPkRatio_BAHAMAS_T8.0_WMAP9_L400N1024.dat");
	sprintf(bary.scenario,"BAHAMAS_T8p0");
	bary.isPkbary = 1 ;
	bary.Nkbins = line_count(bary.FILE_logPkR)-1 ;
	bary.Nabins = 15;		
	double z[15] = {3.0,2.75,2.50,2.25,2.00,1.75,1.50,1.25,1.00,0.75,0.50,0.375,0.25,0.125,0.0};
		
	for (int i=0;i < bary.Nabins;i++) bary.z_bins[i]=z[i];
}

void set_bary_parameters_to_antilles(int id)
{
    if (id < 1 || id > 400) {
        fprintf(stderr, "ANTILLES simulation ID must be 1..400\n");
        exit(EXIT_FAILURE);
    }

    /* Prepare locally; update global state only after validation. */
    barypara next = bary;

    int n = snprintf(next.FILE_logPkR, sizeof next.FILE_logPkR,
        "../cosmolike_core/logPkRatio/logPkRatio_antilles-%d.dat", id);

    if (n < 0 || (size_t)n >= sizeof next.FILE_logPkR) {
        fprintf(stderr, "ANTILLES table path is too long\n");
        exit(EXIT_FAILURE);
    }

    FILE *fp = fopen(next.FILE_logPkR, "r");
    if (fp == NULL) {
        perror(next.FILE_logPkR);
        exit(EXIT_FAILURE);
    }

    char header[1024];
    if (fgets(header, sizeof header, fp) == NULL ||
        strchr(header, '\n') == NULL) {
        fprintf(stderr, "Invalid ANTILLES header: %s\n",
                next.FILE_logPkR);
        fclose(fp);
        exit(EXIT_FAILURE);
    }
    fclose(fp);

    char *token = strtok(header, " \t\r\n");
    if (token == NULL || strcmp(token, "logk") != 0)
        goto bad_header;

    for (int j = 0; j < 15; ++j) {
        token = strtok(NULL, " \t\r\n");
        if (token == NULL || strncmp(token, "z=", 2) != 0)
            goto bad_header;

        char *end;
        double z = strtod(token + 2, &end);
        if (end == token + 2 || *end != '\0' ||
            !isfinite(z) || z < 0.0 ||
            (j > 0 && z >= next.z_bins[j - 1]))
            goto bad_header;

        next.z_bins[j] = z;
    }

    if (strtok(NULL, " \t\r\n") != NULL)
        goto bad_header;

    snprintf(next.scenario, sizeof next.scenario, "antilles-%d", id);
    next.Nkbins = 400;
    next.Nabins = 15;
    next.isPkbary = 1;

    bary = next;
    reset_PkRatio_baryons();
    return;

bad_header:
    fprintf(stderr, "Invalid ANTILLES redshift grid: %s\n",
            next.FILE_logPkR);
    exit(EXIT_FAILURE);
}