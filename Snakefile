rule all:
	input:
		matrix_c = "results/matrix_country_positives.tsv",
		matrix_s = "results/matrix_states_positives.tsv",
		matrix_l = "results/matrix_location_positives.tsv",
		sgtf_c1 = "results/matrix_country_sgtf_growth.tsv",
		sgtf_c2 = "results/matrix_country_sgtf_growth_week.tsv",
		sgtf_s1 = "results/matrix_states_sgtf_growth.tsv",
		sgtf_s2 = "results/matrix_states_sgtf_growth_week.tsv",
		posrate_c1 = "results/matrix_country_positivity_week.tsv",
		posrate_s1 = "results/matrix_states_positivity_week.tsv",
		groups = "results/matrix_agegroups.tsv"

rule arguments:
	params:
		datadir = "data",
		rename_file = "data/rename_columns.xlsx",
		correction_file = "data/fix_values.xlsx",
		cache = "data/combined_testdata.tsv",
		shapefile = "config/bra_admbnda_adm2_ibge_2020.shp",
		coordinates = "config/cache_coordinates.tsv",
		age_groups = "config/demo_bins.txt",
		date_column = "date_testing",
		config_barplot = "/figures/barplot/config_bardetection.tsv",
		config_heatmap1 = "/figures/barplot/config_posrate_weeks_SC2states.tsv",
		config_heatmap2 = "/figures/barplot/config_posrate_weeks_SC2demog.tsv",
		config_lineplot = "/figures/barplot/config_linesgtf.tsv",
		config_map = "/figures/barplot/config_states_sgtf.tsv",
		start_date = "2021-12-01",
		target_week = "2022-08-13"


arguments = rules.arguments.params


rule reformat_hlagyn:
	message:
		"""
		Combine data from HLAGyn
		"""
	input:
		rename = arguments.rename_file,
		correction = "data/fix_values.xlsx",
		cache = arguments.cache
	params:
		datadir = arguments.datadir
	output:
		matrix = "results/combined_testdata_1.tsv",
	shell:
		"""
		python3 scripts/reformat_hlagyn.py \
			--datadir {params.datadir} \
			--rename {input.rename} \
			--correction {input.correction} \
			--cache {input.cache} \
			--output {output.matrix}
		"""


rule reformat_dasa:
	message:
		"""
		Combine data from Dasa
		"""
	input:
		rename = arguments.rename_file,
		correction = "data/fix_values.xlsx",
		cache = rules.reformat_hlagyn.output.matrix
	params:
		datadir = arguments.datadir
	output:
		matrix = "results/combined_testdata_2.tsv",
	shell:
		"""
		python3 scripts/reformat_dasa.py \
			--datadir {params.datadir} \
			--rename {input.rename} \
			--correction {input.correction} \
			--cache {input.cache} \
			--output {output.matrix}
		"""


rule reformat_db:
	message:
		"""
		Combine data from DB Molecular
		"""
	input:
		rename = arguments.rename_file,
		correction = "data/fix_values.xlsx",
		cache = rules.reformat_dasa.output.matrix
	params:
		datadir = arguments.datadir
	output:
		matrix = "results/combined_testdata_3.tsv",
	shell:
		"""
		python3 scripts/reformat_db.py \
			--datadir {params.datadir} \
			--rename {input.rename} \
			--correction {input.correction} \
			--cache {input.cache} \
			--output {output.matrix}
		"""


rule geomatch:
	message:
		"""
		Match location names with geographic shapefile polygons
		"""
	input:
		input_file =  rules.reformat_db.output.matrix,
		coordinates = arguments.coordinates,
		shapefile = arguments.shapefile,
		macros = "config/tabela_municipio_macsaud_estado_regiao.tsv"
	params:
		geo_columns = "state, location",
		add_geo = "country:Brazil",
		lat = "lat",
		long = "long",
		check_match = "ADM2_PT",
		target = "ADM1_PT, ADM1_PCODE, ADM2_PT, ADM2_PCODE",
		target2 = "DS_UF_SIGLA#9, CO_MACSAUD#17, DS_NOMEPAD_macsaud#18, region#7",
		index = "ADM2_PCODE",
		action = "add",
		mode = "columns"
	output:
		matrix = "results/combined_testdata_4.tsv"
	shell:
		"""
		python3 scripts/name2shape.py \
			--input {input.input_file} \
			--shapefile \"{input.shapefile}\" \
			--geo-columns \"{params.geo_columns}\" \
			--add-geo {params.add_geo} \
			--lat {params.lat} \
			--long {params.long} \
			--cache {input.coordinates} \
			--check-match {params.check_match} \
			--target \"{params.target}\" \
			--output {output.matrix}
		
		python3 scripts/reformat_dataframe.py \
			--input1 {output.matrix} \
			--input2 {input.macros} \
			--index {params.index} \
			--action {params.action} \
			--mode {params.mode} \
			--targets "{params.target2}" \
			--output {output.matrix}
		"""


rule demographics:
	message:
		"""
		Aggregate ages per age groups
		"""
	input:
		metadata = rules.geomatch.output.matrix,
		bins = arguments.age_groups,
	params:
		column = "age",
		group = "age_group",
		lowest = "0",
		highest = "200",
		xvar = "age_group",
		format = "integer",
		yvar = "sex geneS_detection epiweek",
		unique_id = "sex",
		filters = "SC2_test_result:Positive, sex:F, sex:M, ~age_group:''",
	output:
		matrix = "results/combined_testdata.tsv",
		age_groups = rules.all.input.groups
	shell:
		"""
		python3 scripts/groupbyrange.py \
			--input {input.metadata} \
			--column {params.column} \
			--bins {input.bins} \
			--group {params.group} \
			--lowest {params.lowest} \
			--highest {params.highest} \
			--output {output.matrix}
			
		python3 scripts/rows2matrix.py \
			--input {output.matrix} \
			--xvar {params.xvar} \
			--format {params.format} \
			--yvar {params.yvar} \
			--filter ""{params.filters}"" \
			--unique-id {params.unique_id} \
			--output {output.age_groups}
		
		"""


rule detection:
	message:
		"""
		Aggregate data related to SGTF results
		"""
	input:
		input_file = "results/combined_testdata.tsv"
	params:
		xvar = arguments.date_column,
		xtype = "time",
		format = "integer",
		
		yvar_country = "country geneS_detection",
		index_country = "country",
		
		yvar_macros = "CO_MACSAUD geneS_detection",
		index_macros = "CO_MACSAUD",
		extra_columns_macros = "DS_UF_SIGLA DS_NOMEPAD_macsaud",

		yvar_states = "ADM1_PCODE geneS_detection",
		index_states = "ADM1_PCODE",
		extra_columns_states = "ADM1_PT DS_UF_SIGLA country",
		
		yvar_location = "ADM2_PCODE geneS_detection",
		index_location = "ADM2_PCODE",
		extra_columns_location = "ADM2_PT state",
		
		filters = "SC2_test_result:Positive, test_kit:thermo",
		start_date = arguments.start_date,
		
		xvar2 = "country",
		extra_cols2 = "ADM1_PT DS_UF_SIGLA",
		time_var = "date_testing",
		start_date2 = "2022-05-15",
		filters_notdetected = "SC2_test_result:Positive, geneS_detection:Not detected, test_kit:thermo",
		filters_detected = "SC2_test_result:Positive, geneS_detection:Detected, test_kit:thermo",
		extra_cols3 = "ADM2_PT state DS_UF_SIGLA lat long",

	output:
		matrix_country = "results/matrix_country_detection.tsv",
		matrix_macros = "results/matrix_macros_detection.tsv",
		matrix_states = "results/matrix_states_detection.tsv",
		matrix_location = "results/matrix_location_detection.tsv",
		pinpoints1 = "figures/map/pinpoints_sgtf.tsv",
		pinpoints2 = "figures/map/pinpoints_sgtp.tsv",
		choropleth1 = "figures/map/choropleth_sgtf.tsv",
		choropleth2 = "figures/map/choropleth_sgtp.tsv"
	shell:
		"""
		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_country} \
			--unique-id {params.index_country} \
			--filter "{params.filters}" \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_country}

		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_macros} \
			--unique-id {params.index_macros} \
			--extra-columns {params.extra_columns_macros} \
			--filter "{params.filters}" \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_macros}
			
		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_states} \
			--unique-id {params.index_states} \
			--extra-columns  {params.extra_columns_states} \
			--filter "{params.filters}" \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_states}
			
		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_location} \
			--unique-id {params.index_location} \
			--extra-columns  {params.extra_columns_location} \
			--filter "{params.filters}" \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_location}
		
		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar2} \
			--format {params.format} \
			--yvar {params.yvar_states} \
			--unique-id {params.index_states} \
			--extra-columns {params.extra_cols2} \
			--filter \"{params.filters_notdetected}\" \
			--time-var {params.time_var} \
			--start-date {params.start_date2} \
			--end-date {arguments.target_week} \
			--output {output.choropleth1}

		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar2} \
			--format {params.format} \
			--yvar {params.yvar_states} \
			--unique-id {params.index_states} \
			--extra-columns {params.extra_cols2} \
			--filter \"{params.filters_detected}\" \
			--time-var {params.time_var} \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.choropleth2}

		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar2} \
			--format {params.format} \
			--yvar {params.yvar_location} \
			--unique-id {params.index_location} \
			--extra-columns {params.extra_cols3} \
			--filter \"{params.filters_notdetected}\" \
			--time-var {params.time_var} \
			--start-date {params.start_date2} \
			--end-date {arguments.target_week} \
			--output {output.pinpoints1}

		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar2} \
			--format {params.format} \
			--yvar {params.yvar_location} \
			--unique-id {params.index_location} \
			--extra-columns {params.extra_cols3} \
			--filter \"{params.filters_detected}\" \
			--time-var {params.time_var} \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.pinpoints2}

		"""


rule test_results:
	message:
		"""
		Aggregate counts of all tests (Postive and Negative)
		"""
	input:
		input_file = "results/combined_testdata.tsv"
	params:
		xvar = arguments.date_column,
		xtype = "time",
		format = "integer",
		
		yvar_country = "country SC2_test_result",
		index_country = "country",
		
		yvar_states = "DS_UF_SIGLA SC2_test_result",
		index_states = "DS_UF_SIGLA",
		extra_columns_states = "ADM1_PT country",
		
		yvar_location = "ADM2_PCODE SC2_test_result",
		index_location = "ADM2_PCODE",
		extra_columns_location = "ADM2_PT state",
		
		start_date = arguments.start_date,
		filters = "~SC2_test_result:Não detectado , ~SC2_test_result:Inconclusivo"
	output:
		matrix_country = "results/matrix_country_posneg_all.tsv",
		matrix_states = "results/matrix_states_posneg_all.tsv",
		matrix_location = "results/matrix_location_posneg_all.tsv",
	shell:
		"""
		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_country} \
			--unique-id {params.index_country} \
			--filter \"{params.filters}\" \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_country}

		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_states} \
			--unique-id {params.index_states} \
			--filter \"{params.filters}\" \
			--extra-columns  {params.extra_columns_states} \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_states}
			
		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_location} \
			--unique-id {params.index_location} \
			--filter \"{params.filters}\" \
			--extra-columns  {params.extra_columns_location} \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_location}
		"""


rule proportions:
	message:
		"""
		Proportions of positives and SGTF
		"""
	input:
		file = "results/combined_testdata.tsv",
	params:
		xvar = "DS_UF_SIGLA",
		format = "integer",
		filters = "test_kit:thermo",
		
		yvar1 = "epiweek geneS_detection",
		index1 = "epiweek",
		
		yvar2 = "epiweek SC2_test_result",
		index2 = "epiweek",
	output:
		file1 = "results/matrix_states_proportions_sgtf_weeks.tsv",
		file2 = "results/matrix_states_proportions_posneg_weeks.tsv"
	shell:
		"""
		python3 scripts/rows2matrix.py \
			--input {input.file} \
			--xvar {params.xvar} \
			--format {params.format} \
			--yvar {params.yvar1} \
			--unique-id {params.index1} \
			--filter "{params.filters}" \
			--output {output.file1}

		python3 scripts/rows2matrix.py \
			--input {input.file} \
			--xvar {params.xvar} \
			--format {params.format} \
			--yvar {params.yvar2} \
			--unique-id {params.index2} \
			--filter "{params.filters}" \
			--output {output.file2}
		"""


rule aggregate:
	message:
		"""
		Aggregate data by week
		"""
	input:
		input_c1 = "results/matrix_country_detection.tsv",
		input_c2 = "results/matrix_country_positives.tsv",
		input_c3 = "results/matrix_country_posneg_all.tsv",
		input_s1 = "results/matrix_states_detection.tsv",
		input_s2 = "results/matrix_states_positives.tsv",
		input_s3 = "results/matrix_states_posneg_all.tsv",
	params:
		unit = "week",
		format = "integer",
		week_end = "end",
	output:
		matrix_c1 = "results/matrix_country_detection_week.tsv",
		matrix_c2 = "results/matrix_country_positives_week.tsv",
		matrix_c3 = "results/matrix_country_posneg_all_week.tsv",
		matrix_s1 = "results/matrix_states_detection_week.tsv",
		matrix_s2 = "results/matrix_states_positives_week.tsv",
		matrix_s3 = "results/matrix_states_posneg_all_week.tsv",

	shell:
		"""
		python3 scripts/aggregator.py \
			--input {input.input_c1} \
			--unit {params.unit} \
			--weekasdate {params.week_end} \
			--format {params.format} \
			--output {output.matrix_c1}
		
		python3 scripts/aggregator.py \
			--input {input.input_c2} \
			--unit {params.unit} \
			--weekasdate {params.week_end} \
			--format {params.format} \
			--output {output.matrix_c2}
		
		python3 scripts/aggregator.py \
			--input {input.input_c3} \
			--unit {params.unit} \
			--weekasdate {params.week_end} \
			--format {params.format} \
			--output {output.matrix_c3}

		python3 scripts/aggregator.py \
			--input {input.input_s1} \
			--unit {params.unit} \
			--weekasdate {params.week_end} \
			--format {params.format} \
			--output {output.matrix_s1}
		
		python3 scripts/aggregator.py \
			--input {input.input_s2} \
			--unit {params.unit} \
			--weekasdate {params.week_end} \
			--format {params.format} \
			--output {output.matrix_s2}

		python3 scripts/aggregator.py \
			--input {input.input_s3} \
			--unit {params.unit} \
			--weekasdate {params.week_end} \
			--format {params.format} \
			--output {output.matrix_s3}
		
		"""



rule sgtf_percent:
	message:
		"""
		Percentage of STGF cases
		"""
	input:
		input_c1 = "results/matrix_country_detection.tsv",
		input_c2 = "results/matrix_country_positives.tsv",
		input_c3 = "results/matrix_country_detection_week.tsv",
		input_c4 = "results/matrix_country_positives_week.tsv",
		input_s1 = "results/matrix_states_detection.tsv",
		input_s2 = "results/matrix_states_positives.tsv",
		input_s3 = "results/matrix_states_detection_week.tsv",
		input_s4 = "results/matrix_states_positives_week.tsv",
	params:
		index_c1 = "country geneS_detection",
		index_c2 = "country",
		index_s1 = "ADM1_PCODE geneS_detection",
		index_s2 = "ADM1_PCODE"
	output:
		matrix_c1 = rules.all.input.sgtf_c1,
		matrix_c2 = rules.all.input.sgtf_c2,
		matrix_s1 = rules.all.input.sgtf_s1,
		matrix_s2 = rules.all.input.sgtf_s2,
	shell:
		"""
		python3 scripts/matrix_operations.py \
			--input1 {input.input_c1} \
			--input2 {input.input_c2} \
			--index1 {params.index_c1} \
			--index2 {params.index_c2} \
			--output {output.matrix_c1}

		python3 scripts/matrix_operations.py \
			--input1 {input.input_c3} \
			--input2 {input.input_c4} \
			--index1 {params.index_c1} \
			--index2 {params.index_c2} \
			--output {output.matrix_c2}

		python3 scripts/matrix_operations.py \
			--input1 {input.input_s1} \
			--input2 {input.input_s2} \
			--index1 {params.index_s1} \
			--index2 {params.index_s2} \
			--output {output.matrix_s1}

		python3 scripts/matrix_operations.py \
			--input1 {input.input_s3} \
			--input2 {input.input_s4} \
			--index1 {params.index_s1} \
			--index2 {params.index_s2} \
			--output {output.matrix_s2}
		
		"""


rule all_positives:
	message:
		"""
		Aggregate counts of all positive SGTF tests
		"""
	input:
		input_file = "results/combined_testdata.tsv"
	params:
		xvar = arguments.date_column,
		xtype = "time",
		format = "integer",

		yvar_country = "country",
		index_country = "country",
		
		yvar_states = "ADM1_PCODE",
		index_states = "ADM1_PCODE",
		extra_columns_states = "ADM1_PT country",
		
		yvar_location = "ADM2_PCODE",
		index_location = "ADM2_PCODE",
		extra_columns_location = "ADM2_PT state",
		
		filters = "SC2_test_result:Positive, test_kit:thermo",
		start_date = arguments.start_date,
	output:
		matrix_country = rules.all.input.matrix_c,
		matrix_states = rules.all.input.matrix_s,
		matrix_location = rules.all.input.matrix_l,
	shell:
		"""
		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_country} \
			--unique-id {params.index_country} \
			--filter "{params.filters}" \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_country}

		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_states} \
			--unique-id {params.index_states} \
			--extra-columns  {params.extra_columns_states} \
			--filter "{params.filters}" \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_states}
			
		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_location} \
			--unique-id {params.index_location} \
			--extra-columns  {params.extra_columns_location} \
			--filter "{params.filters}" \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_location}
		"""


rule total_tests:
	message:
		"""
		Aggregate counts of total SGTF tests
		"""
	input:
		input_file = "results/combined_testdata.tsv"
	params:
		xvar = arguments.date_column,
		xtype = "time",
		format = "integer",
		filters = "test_kit:thermo",

		yvar_country = "country",
		index_country = "country",
		
		yvar_states = "ADM1_PCODE",
		index_states = "ADM1_PCODE",
		extra_columns_states = "ADM1_PT country",
		
		yvar_location = "ADM2_PCODE",
		index_location = "ADM2_PCODE",
		extra_columns_location = "ADM2_PT state",
		
		start_date = arguments.start_date,
	output:
		matrix_country = "results/matrix_country_total.tsv",
		matrix_states = "results/matrix_states_total.tsv",
		matrix_location = "results/matrix_location_total.tsv",
		matrix_country2 = "results/matrix_country_total_all.tsv",
		matrix_states2 = "results/matrix_states_total_all.tsv",
	shell:
		"""
		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_country} \
			--unique-id {params.index_country} \
			--filter "{params.filters}" \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_country}

		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_states} \
			--unique-id {params.index_states} \
			--extra-columns  {params.extra_columns_states} \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--filter "{params.filters}" \
			--output {output.matrix_states}
			
		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_location} \
			--unique-id {params.index_location} \
			--extra-columns  {params.extra_columns_location} \
			--filter "{params.filters}" \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_location}

		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_country} \
			--unique-id {params.index_country} \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_country2}

		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--xtype {params.xtype} \
			--format {params.format} \
			--yvar {params.yvar_states} \
			--unique-id {params.index_states} \
			--extra-columns  {params.extra_columns_states} \
			--start-date {params.start_date} \
			--end-date {arguments.target_week} \
			--output {output.matrix_states2}

		"""



rule posrate_agegroup:
	message:
		"""
		Positive rate by age group
		"""
	input:
		input_file = "results/combined_testdata.tsv"
	params:
		format = "integer",
		xvar = "epiweek",
		yvar = "age_group SC2_test_result",
		unique_id = "age_group",
		extra = "country",
		min_denominator = 50,
		filters = "~SC2_test_result:Não detectado , ~SC2_test_result:Inconclusivo"
	output:
		week_matrix = "results/matrix_agegroups_weeks_SC2_posneg.tsv",
		alltests = "results/matrix_agegroups_weeks_SC2_alltests.tsv",
		posrate = "results/matrix_agegroups_weeks_SC2_posrate.tsv",
	shell:
		"""
		python3 scripts/rows2matrix.py \
			--input {input.input_file} \
			--xvar {params.xvar} \
			--format {params.format} \
			--yvar {params.yvar} \
			--extra-columns {params.extra} \
			--unique-id {params.unique_id} \
			--filter "{params.filters}" \
			--output {output.week_matrix}

		python3 scripts/collapser.py \
			--input {output.week_matrix} \
			--index {params.unique_id} \
			--unique-id {params.unique_id} \
			--extra-columns {params.extra} \
			--output {output.alltests} \

		python3 scripts/matrix_operations.py \
			--input1 {output.week_matrix} \
			--input2 {output.alltests} \
			--index1 {params.yvar} \
			--index2 {params.unique_id} \
			--min-denominator {params.min_denominator} \
			--output {output.posrate}
		"""


rule positivity:
	message:
		"""
		Percentage of positive tests
		"""
	input:
		input_c1 = "results/matrix_country_posneg_all_week.tsv",
		input_s1 = "results/matrix_states_posneg_all_week.tsv",
	params:
		index_c1 = "country SC2_test_result",
		index_c2 = "country",
		index_s1 = "DS_UF_SIGLA SC2_test_result",
		index_s2 = "DS_UF_SIGLA",
		format = "integer",
		ignore_s = "ADM1_PT country",
		min_denom = "50"
	output:
		matrix_cd = "results/matrix_country_alltests_week.tsv",
		matrix_sd = "results/matrix_states_alltests_week.tsv",
		matrix_c1 = rules.all.input.posrate_c1,
		matrix_s1 = rules.all.input.posrate_s1,
	shell:
		"""
		python3 scripts/collapser.py \
			--input {input.input_c1} \
			--index {params.index_c2} \
			--unique-id {params.index_c2} \
			--format {params.format} \
			--output {output.matrix_cd}

		python3 scripts/collapser.py \
			--input {input.input_s1} \
			--index {params.index_s2} \
			--unique-id {params.index_s2} \
			--ignore {params.ignore_s} \
			--format {params.format} \
			--output {output.matrix_sd}

		python3 scripts/matrix_operations.py \
			--input1 {input.input_c1} \
			--input2 {output.matrix_cd} \
			--index1 {params.index_c1} \
			--index2 {params.index_c2} \
			--min-denominator {params.min_denom} \
			--output {output.matrix_c1}

		python3 scripts/matrix_operations.py \
			--input1 {input.input_s1} \
			--input2 {output.matrix_sd} \
			--index1 {params.index_s1} \
			--index2 {params.index_s2} \
			--min-denominator {params.min_denom} \
			--output {output.matrix_s1}
		
		"""






rule dataviz:
	message:
		"""
		Create data visualizations
		"""
	input:
		barplot = arguments.config_barplot,
		heatmap1 = arguments.config_heatmap1,
		heatmap2 = arguments.config_heatmap2,
		lineplot = arguments.config_lineplot,
		map = arguments.config_map,
	shell:
		"""
		python3 scripts/pandas_multibar.py \
			--config {input.barplot}

		python3 scripts/pandas_heatmap.py \
			--config {input.heatmap1}

		python3 scripts/pandas_heatmap.py \
			--config {input.heatmap2}

		python3 scripts/pandas_lineplot.py \
			--config {input.lineplot}

		python3 scripts/duomap.py \
			--config {input.map}
		"""


rule stats:
	message:
		"""
		Generate file with summary statistics
		"""
	input:
		input_file = "results/combined_testdata.tsv"
	params:
		week = arguments.target_week
	output:
		counts = "results/counts.txt"
	shell:
		"""
		echo "Test positivity" > {output.counts}
		cut -d$'\t' -f 10 {input.input_file} | grep -v SC2_test_result | sort | uniq -c >> {output.counts}

		echo "\nNumber of states included in the surveillance" >> {output.counts}
		cut -d$'\t' -f 24 {input.input_file} | sort | uniq -c | wc -l >> {output.counts}

		echo "\nNumber of macroregions of health included in the surveillance" >> {output.counts}
		cut -d$'\t' -f 24,26 {input.input_file} | sort | uniq -c | wc -l >> {output.counts}

		echo "\nNumber of municipalities included in the surveillance" >> {output.counts}
		cut -d$'\t' -f 24,22 {input.input_file} | sort | uniq -c | wc -l >> {output.counts}

		echo "\nPrevalence SGTF Brazil" >> {output.counts}
		cut -d$'\t' -f 10,16,13 {input.input_file} | sort | grep Positive | grep {params.week} | uniq -c >> {output.counts} 

		echo "\nPositive and Negative tests in the last EW in Brazil" >> {output.counts}
		cut -d$'\t' -f 10,13 {input.input_file} | sort | grep {params.week} | uniq -c >> {output.counts} 

		echo "\nStates with cases of Omicron" >> {output.counts}
		cut -d$'\t' -f 10,16,24 {input.input_file} | sort | grep -v Negative | grep detected | uniq >> {output.counts} 

		echo "\nMacroregions with cases of Omicron" >> {output.counts}
		cut -d$'\t' -f 10,16,24,26 {input.input_file} | sort | grep -v Negative | grep detected | uniq -c | wc -l >> {output.counts} 

		echo "\nNumber of municipalities with cases of Omicron" >> {output.counts}
		cut -d$'\t' -f 10,16,24,22 {input.input_file} | sort | grep -v Negative | grep detected | grep -v SC2_test_result | uniq | wc -l >> {output.counts}
		"""


rule copy_files:
	message:
		"""
		Copy files for plotting
		"""
	shell:
		"""
		cp "results/matrix_country_detection_week.tsv" ./figures/barplot
		cp "results/matrix_country_sgtf_growth_week.tsv" ./figures/lineplot

		"""


rule remove_figs:
	message: "Removing figures"
	shell:
		"""
		rm figures/*/*.pdf
		rm figures/*/matrix*
		rm figures/*/pinpoints*
		rm figures/*/choropleth*
		"""



#
#
#rule xxx:
#	message:
#		"""
#		
#		"""
#	input:
#		metadata = arguments.
#	params:
#		index = arguments.,
#		date = arguments.
#	output:
#		matrix = "results/"
#	shell:
#		"""
#		python3 scripts/ \
#			--metadata {input.} \
#			--index-column {params.} \
#			--extra-columns {params.} \
#			--date-column {params.} \
#			--output {output.}
#		"""

rule clean:
	message: "Removing directories: {params}"
	params:
		"results"
	shell:
		"""
		rm -rfv {params}
		"""
