Search.setIndex({docnames:["changelog","cli","contributing","contributors","format","genindex","index","introduction","license","modules","robocrys","robocrys.condense","robocrys.describe","robocrys.featurize"],envversion:{"sphinx.domains.c":2,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":3,"sphinx.domains.index":1,"sphinx.domains.javascript":2,"sphinx.domains.math":2,"sphinx.domains.python":2,"sphinx.domains.rst":2,"sphinx.domains.std":2,"sphinx.ext.intersphinx":1,"sphinx.ext.todo":2,"sphinx.ext.viewcode":1,sphinx:56},filenames:["changelog.rst","cli.rst","contributing.rst","contributors.rst","format.rst","genindex.rst","index.rst","introduction.rst","license.rst","modules.rst","robocrys.rst","robocrys.condense.rst","robocrys.describe.rst","robocrys.featurize.rst"],objects:{"":{robocrys:[10,0,0,"-"]},"robocrys.StructureCondenser":{condense_structure:[10,2,1,""]},"robocrys.StructureDescriber":{describe:[10,2,1,""],get_all_component_descriptions:[10,2,1,""],get_bond_length_description:[10,2,1,""],get_component_description:[10,2,1,""],get_component_makeup_summary:[10,2,1,""],get_mineral_description:[10,2,1,""],get_octahedral_tilt_description:[10,2,1,""],get_site_description:[10,2,1,""]},"robocrys.adapter":{BaseAdapter:[10,1,1,""]},"robocrys.adapter.BaseAdapter":{angles:[10,2,1,""],component_makeup:[10,2,1,""],components:[10,2,1,""],crystal_system:[10,2,1,""],dimensionality:[10,2,1,""],distances:[10,2,1,""],elements:[10,3,1,""],formula:[10,2,1,""],get_angle_details:[10,2,1,""],get_distance_details:[10,2,1,""],is_vdw_heterostructure:[10,2,1,""],mineral:[10,2,1,""],sites:[10,2,1,""],spg_symbol:[10,2,1,""]},"robocrys.cli":{main:[10,4,1,""],robocrystallographer:[10,4,1,""]},"robocrys.condense":{component:[11,0,0,"-"],condenser:[11,0,0,"-"],fingerprint:[11,0,0,"-"],mineral:[11,0,0,"-"],molecule:[11,0,0,"-"],site:[11,0,0,"-"]},"robocrys.condense.component":{components_are_isomorphic:[11,4,1,""],components_are_vdw_heterostructure:[11,4,1,""],filter_molecular_components:[11,4,1,""],get_component_formula:[11,4,1,""],get_component_formula_and_factor:[11,4,1,""],get_formula_from_components:[11,4,1,""],get_formula_inequiv_components:[11,4,1,""],get_reconstructed_structure:[11,4,1,""],get_structure_inequiv_components:[11,4,1,""],get_sym_inequiv_components:[11,4,1,""],get_vdw_heterostructure_information:[11,4,1,""]},"robocrys.condense.condenser":{StructureCondenser:[11,1,1,""]},"robocrys.condense.condenser.StructureCondenser":{condense_structure:[11,2,1,""]},"robocrys.condense.fingerprint":{get_fingerprint_distance:[11,4,1,""],get_site_fingerprints:[11,4,1,""],get_structure_fingerprint:[11,4,1,""]},"robocrys.condense.mineral":{MineralMatcher:[11,1,1,""]},"robocrys.condense.mineral.MineralMatcher":{get_aflow_matches:[11,2,1,""],get_best_mineral_name:[11,2,1,""],get_fingerprint_matches:[11,2,1,""]},"robocrys.condense.molecule":{MoleculeNamer:[11,1,1,""]},"robocrys.condense.molecule.MoleculeNamer":{get_name_from_molecule_graph:[11,2,1,""],get_name_from_pubchem:[11,2,1,""],molecule_graph_to_smiles:[11,2,1,""],name_sources:[11,3,1,""]},"robocrys.condense.site":{SiteAnalyzer:[11,1,1,""],geometries_match:[11,4,1,""],nn_summaries_match:[11,4,1,""],nnn_summaries_match:[11,4,1,""]},"robocrys.condense.site.SiteAnalyzer":{equivalent_sites:[11,3,1,""],get_all_bond_distance_summaries:[11,2,1,""],get_all_connectivity_angle_summaries:[11,2,1,""],get_all_nnn_distance_summaries:[11,2,1,""],get_all_site_summaries:[11,2,1,""],get_bond_distance_summary:[11,2,1,""],get_connectivity_angle_summary:[11,2,1,""],get_inequivalent_site_indices:[11,2,1,""],get_nearest_neighbors:[11,2,1,""],get_next_nearest_neighbors:[11,2,1,""],get_nnn_distance_summary:[11,2,1,""],get_site_geometry:[11,2,1,""],get_site_summary:[11,2,1,""],symmetry_labels:[11,3,1,""]},"robocrys.describe":{adapter:[12,0,0,"-"],describer:[12,0,0,"-"]},"robocrys.describe.adapter":{ComponentDetails:[12,1,1,""],ComponentGroup:[12,1,1,""],DescriptionAdapter:[12,1,1,""],NeighborSiteDetails:[12,1,1,""],NextNeighborSiteDetails:[12,1,1,""],SiteGroup:[12,1,1,""]},"robocrys.describe.adapter.ComponentDetails":{count:[12,3,1,""],dimensionality:[12,3,1,""],formula:[12,3,1,""],index:[12,3,1,""],molecule_name:[12,3,1,""],nsites:[12,3,1,""],orientation:[12,3,1,""]},"robocrys.describe.adapter.ComponentGroup":{components:[12,3,1,""],count:[12,3,1,""],dimensionality:[12,3,1,""],formula:[12,3,1,""],molecule_name:[12,3,1,""],nsites:[12,3,1,""]},"robocrys.describe.adapter.DescriptionAdapter":{get_component_details:[12,2,1,""],get_component_groups:[12,2,1,""],get_component_site_groups:[12,2,1,""],get_nearest_neighbor_details:[12,2,1,""],get_next_nearest_neighbor_details:[12,2,1,""],get_sym_label:[12,2,1,""],sym_labels:[12,3,1,""],use_iupac_ordering:[12,3,1,""]},"robocrys.describe.adapter.NeighborSiteDetails":{count:[12,3,1,""],element:[12,3,1,""],sites:[12,3,1,""],sym_label:[12,3,1,""]},"robocrys.describe.adapter.NextNeighborSiteDetails":{connectivity:[12,3,1,""],count:[12,3,1,""],element:[12,3,1,""],geometry:[12,3,1,""],poly_formula:[12,3,1,""],sites:[12,3,1,""],sym_label:[12,3,1,""]},"robocrys.describe.adapter.SiteGroup":{count:[12,3,1,""],element:[12,3,1,""],sites:[12,3,1,""]},"robocrys.describe.describer":{StructureDescriber:[12,1,1,""],get_mineral_name:[12,4,1,""]},"robocrys.describe.describer.StructureDescriber":{describe:[12,2,1,""],get_all_component_descriptions:[12,2,1,""],get_bond_length_description:[12,2,1,""],get_component_description:[12,2,1,""],get_component_makeup_summary:[12,2,1,""],get_mineral_description:[12,2,1,""],get_octahedral_tilt_description:[12,2,1,""],get_site_description:[12,2,1,""]},"robocrys.featurize":{adapter:[13,0,0,"-"],featurizer:[13,0,0,"-"]},"robocrys.featurize.adapter":{FeaturizerAdapter:[13,1,1,""]},"robocrys.featurize.adapter.FeaturizerAdapter":{all_bond_lengths:[13,2,1,""],average_anion_coordination_number:[13,2,1,""],average_cation_coordination_number:[13,2,1,""],average_coordination_number:[13,2,1,""],average_corner_sharing_octahedral_tilt_angle:[13,2,1,""],component_dimensionalities:[13,2,1,""],contains_connected_geometry:[13,2,1,""],contains_corner_sharing_polyhedra:[13,2,1,""],contains_edge_sharing_polyhedra:[13,2,1,""],contains_face_sharing_polyhedra:[13,2,1,""],contains_geometry_type:[13,2,1,""],contains_molecule:[13,2,1,""],contains_named_molecule:[13,2,1,""],contains_polyhedra:[13,2,1,""],frac_site_geometry:[13,2,1,""],frac_sites_n_coordinate:[13,2,1,""],frac_sites_polyhedra:[13,2,1,""],is_dimensionality:[13,2,1,""],is_intercalated:[13,2,1,""],is_interpenetrated:[13,2,1,""]},"robocrys.featurize.featurizer":{RobocrysFeaturizer:[13,1,1,""]},"robocrys.featurize.featurizer.RobocrysFeaturizer":{citations:[13,2,1,""],feature_labels:[13,2,1,""],featurize:[13,2,1,""],implementors:[13,2,1,""]},"robocrys.util":{common_formulas:[10,3,1,""],connected_geometries:[10,3,1,""],defaultdict_to_dict:[10,4,1,""],dimensionality_to_shape:[10,3,1,""],geometry_to_polyhedra:[10,3,1,""],get_el:[10,4,1,""],get_formatted_el:[10,4,1,""],htmlify_spacegroup:[10,4,1,""],load_condensed_structure_json:[10,4,1,""],superscript_number:[10,4,1,""],unicodeify_spacegroup:[10,4,1,""]},robocrys:{StructureCondenser:[10,1,1,""],StructureDescriber:[10,1,1,""],adapter:[10,0,0,"-"],cli:[10,0,0,"-"],condense:[11,0,0,"-"],describe:[12,0,0,"-"],featurize:[13,0,0,"-"],util:[10,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","class","Python class"],"2":["py","method","Python method"],"3":["py","attribute","Python attribute"],"4":["py","function","Python function"]},objtypes:{"0":"py:module","1":"py:class","2":"py:method","3":"py:attribute","4":"py:function"},terms:{"001":[1,11],"017":11,"100":2,"1016":11,"135":4,"136":11,"1557":[6,7],"2005":[10,11,12],"2015":8,"2017":11,"2018":1,"2019":[6,7],"42_1m":10,"5544":11,"856":[1,6,7],"874":[6,7],"881":[6,7],"case":[10,11],"class":[6,7,10,11,12,13],"default":[1,10,11],"float":[10,11,13],"function":[8,10,11],"import":[6,7,10],"int":[10,11,12,13],"new":2,"public":[6,7],"return":[10,11,12,13],"sno\u2082":10,"static":11,"super":10,"throw":[10,12],"true":[1,4,10,11,12],AND:8,ARE:8,BUT:8,FOR:8,For:[1,2,6,7,10,11],Has:11,MRS:[6,7],NOT:8,One:[10,13],PRs:2,SUCH:8,THE:8,THe:12,The:[1,2,6,7,8,10,11,12,13],Their:10,These:1,USE:8,Use:[2,10],Uses:11,Will:[10,12],about:11,abov:8,absolut:2,accept:2,access:[1,10],accord:[10,11,12],account:[11,13],across:13,actanid:[10,11,12],actual:2,adapt:[2,9],added:[10,12],addit:11,advis:8,aflow:11,aflow_prototyp:11,aflowprototypematch:11,against:11,agreement:8,ajain:13,alex:[1,3],algorithm:11,alia:12,all:[2,6,7,8,10,11,12,13],all_bond_length:13,allow:[1,6,7],along:11,alreadi:[10,12],also:[1,6,7],altern:[1,6,7],alwai:[1,10,11],analys:[6,7],analysi:[6,7,10,11],analyz:11,angl:[1,4,6,7,10,11,12,13],angle_1:11,angle_2:11,ani:[2,8,10,11,12,13],anion:[1,13],anubhav:[3,13],anyon:8,api:[1,6,7],api_kei:[1,6,7],appear:11,appli:10,applic:12,approach:11,approv:8,arbitrari:10,area:2,argument:[10,13],aris:8,as_dict:11,atom:[2,6,7,10],attribut:[12,13],author:[1,13],autom:[6,7],automat:[1,2],avail:8,averag:13,average_anion_coordination_numb:13,average_cation_coordination_numb:13,average_coordination_numb:13,average_corner_sharing_octahedral_tilt_angl:13,avoid:2,axi:11,background:2,base:[10,11,12,13],baseadapt:[10,12,13],basefeatur:13,basic:2,becom:2,been:[10,11,12],begin:11,being:2,below:[1,4,6,7,8],beneath:11,berkelei:[3,8],best:[1,2,6,7,11],between:[10,11,12],bibtex:13,binari:8,bond:[1,6,7,10,11,12],bond_angle_tol:11,bond_dist_tol:11,bond_length_decimal_plac:[10,12],bonded_structur:11,bool:[10,11,12,13],both:[1,6,7,10,11],branch:2,briandk:2,broken:11,bsd:[6,7,8],bug:8,busi:8,calcul:[10,11],california:8,can:[1,2,6,7,10,11,12],career:3,cation:[1,13],caus:8,cell:[1,10,11],centr:11,chang:[2,6,7,10],changelog:[6,7],charact:10,check:11,chemic:11,chemistri:[10,11,12],choos:8,cif:[1,6,7],citat:13,classifi:13,cli:9,close:11,code:8,codebas:2,combin:11,command:[6,7,10],commatsci:11,comment:2,commit:2,common:[1,6,7,10,11],common_formula:10,commun:[6,7],compar:[2,11],compon:[1,4,6,7,9,10,12,13],component_a:11,component_b:11,component_dimension:13,component_index:[10,12],component_makeup:[4,10,12],componentdetail:12,componentgroup:12,components_are_isomorph:11,components_are_vdw_heterostructur:11,composit:[10,11],composition:11,compound:11,comput:[8,11],con:11,conda:[6,7],condens:[6,7,9,10,12,13],condense_structur:[6,7,10,11,12,13],condensed_structur:[6,7,10,12,13],condenser_kwarg:[10,13],condit:8,confirm:2,connect:[1,6,7,10,11,12,13],connected_geometri:[10,11,12],connectivity_a:11,connectivity_b:11,consequenti:8,consid:[6,7,10,11],consist:11,constrain:11,contain:[1,10,11,12,13],contains_connected_geometri:13,contains_corner_sharing_polyhedra:13,contains_edge_sharing_polyhedra:13,contains_face_sharing_polyhedra:13,contains_geometry_typ:13,contains_molecul:13,contains_named_molecul:13,contains_polyhedra:13,content:9,contract:8,contributor:[2,6,7,8],control:1,convent:[1,10,11],convert:[1,10,11,12],coordin:[4,11,13],copyright:8,core:[6,7,10],corner:[4,6,7,10,11,12,13],correspond:11,count:[11,12],creat:2,crystal:[1,6,7,10],crystal_system:[4,10],crystallograph:[1,6,7,11],crystalnn:[10,11],crystalnnfingerprint_op:11,current:[2,10,12],curtarolo:11,customis:[6,7],cutoff:11,damag:8,data:[8,10,11,12,13],databas:[6,7,10,11],date:[10,11],decemb:1,decim:1,decreas:11,dedic:2,defaultdict:10,defaultdict_to_dict:10,defin:[10,11,12],depart:3,depend:11,dept:8,der:11,deriv:[8,10,12],describ:[6,7,9,10],describe_bond_length:[10,12],describe_compon:[10,12],describe_component_makeup:[10,12],describe_miner:[10,12],describe_oxidation_st:[10,12],describe_symmetry_label:[10,12],describer_kwarg:10,descript:[1,6,7,10,12],descriptionadapt:12,detail:[1,6,7,10,11],determin:[1,6,7,10,11],develop:[2,3,6,7,8],dict:[10,11,12,13],dictionari:[10,11,12,13],didn:2,differ:[10,11,12],dimension:[4,6,7,10,11,12,13],dimensionality_to_shap:10,direct:[8,11],directli:[1,6,7,8],disabl:1,disclaim:8,dist:11,dist_1:11,dist_2:11,dist_3:11,distanc:[4,10,11,12],distance_1:11,distance_2:11,distinct:11,distort:[1,4,10,11,12,13],distorted_tol:[10,12,13],distribut:8,doc:1,docs_rst:[10,11],document:[2,6,7,8,10,11],doe:13,doi:[6,7,11],don:[1,10],draft:2,duplic:11,each:[1,10,11,12,13],earli:[3,6,7],easi:2,edg:[4,6,7,10,11,13],effect:[10,11,12],either:[8,13],electron:[1,6,7,11],electroneg:[10,11,12],element:[1,2,4,10,11,12,13],elong:11,els:[11,12],email:13,empti:[10,12],encount:2,endors:8,energi:[3,8],enhanc:8,ensur:2,environ:[6,7],equal:[10,12],equival:[1,6,7,10,11],equivalent_sit:11,error:[2,10,12],etc:[10,11],euclidean:11,even:[1,8],event:8,exampl:[1,2,4,10,11],excel:2,except:[10,11,12],exclus:8,exemplari:8,expect:2,express:8,extend:[6,7],extract:11,face:[4,10,11,13],facebook:2,facilit:[10,12,13],factor:11,fall:[2,11],fals:[1,4,10,11,12],featur:[2,8,9,10,11],feature_label:13,featurizeradapt:13,fetch:[6,7],few:2,field:[11,12],file:[1,6,7,10],filenam:[1,10],filter:11,filter_molecular_compon:11,find:2,fingerprint:[9,10],fingerprint_distance_cutoff:11,fingerprint_tol:11,finish:2,first:[1,11],fit:8,fix:[2,8,10,11],flow:2,fmt:[1,10,12],follow:[2,6,7,8,10,11,12],followup:2,forg:[6,7],fork:2,form:[6,7,8,11],format:[1,10,11,12,13],formula:[1,4,10,11,12],formuula:11,forum:2,found:[6,7,11],frac_site_geometri:13,frac_sites_n_coordin:13,frac_sites_polyhedra:13,fraction:[11,13],framework:2,free:8,from:[1,2,3,6,7,8,10,11,12,13],from_fil:[6,7],from_sit:[10,11,12],full:[1,2,6,7],fund:3,further:2,ga2as2:11,gaa:11,ganos:[1,3,6,7],gener:[1,6,7,10,11,12,13],geom:11,geometr:11,geometri:[1,4,6,7,10,11,12,13],geometries_match:11,geometry_a:11,geometry_b:11,geometry_to_polyhedra:10,geometry_typ:11,get:[10,11,12],get_aflow_match:11,get_all_bond_distance_summari:11,get_all_component_descript:[10,12],get_all_connectivity_angle_summari:11,get_all_nnn_distance_summari:11,get_all_site_summari:11,get_angle_detail:10,get_best_mineral_nam:11,get_bond_distance_summari:11,get_bond_length_descript:[10,12],get_component_descript:[10,12],get_component_detail:12,get_component_formula:11,get_component_formula_and_factor:11,get_component_group:12,get_component_makeup_summari:[10,12],get_component_site_group:12,get_connectivity_angle_summari:11,get_distance_detail:10,get_el:10,get_fingerprint_dist:11,get_fingerprint_match:11,get_formatted_el:10,get_formula_from_compon:11,get_formula_inequiv_compon:11,get_inequivalent_site_indic:11,get_mineral_descript:[10,12],get_mineral_nam:12,get_name_from_molecule_graph:11,get_name_from_pubchem:11,get_nearest_neighbor:11,get_nearest_neighbor_detail:12,get_next_nearest_neighbor:11,get_next_nearest_neighbor_detail:12,get_nnn_distance_summari:11,get_octahedral_tilt_descript:[10,12],get_reconstructed_structur:11,get_site_descript:[10,12],get_site_fingerprint:11,get_site_geometri:11,get_site_summari:11,get_structure_by_material_id:[6,7],get_structure_compon:11,get_structure_fingerprint:11,get_structure_inequiv_compon:11,get_sym_inequiv_compon:11,get_sym_inequiv_equival:11,get_sym_label:12,get_vdw_heterostructure_inform:11,give:11,given:[4,6,7,10,11],good:[2,8,10,12],googl:2,gov:[10,11,13],gradual:11,grant:8,graph:11,group:[2,3,6,7,10,11,12],group_by_el:12,guid:2,guidelin:[2,6,7],gzip:1,hack:3,hackingmateri:[10,11],handl:[11,12],hanson:11,happen:2,hart:11,has:[1,11],have:[2,6,7,10,11,12],here:[2,6,7,11],herebi:8,heteostructur:11,heterostructur:[6,7,10,11],hexagon:4,hick:11,hitchhik:2,holder:8,host:2,how:[1,10,11],howev:[1,8],html:[1,10,11],htmlify_spacegroup:10,http:[1,6,7,10,11],hydrogen:[10,11,12],ideal:13,identifi:[6,7],ids:[1,6,7],ignor:11,implement:[2,10,11,12,13],implementor:13,impli:8,impos:8,improv:2,inc_graph:11,inc_inequivalent_site_index:11,inc_intercal:11,inc_molecule_graph:11,inc_ordered_compon:11,inc_orient:11,inc_site_id:11,incident:8,includ:[1,2,6,7,8,10,11],incorpor:[2,8],index:[10,11,12],indic:[10,11,12],indirect:8,individu:11,inequiv_index:11,inequival:[1,10,11,12],inform:[1,2,6,7,10,11],init:11,initial_angle_tol:11,initial_ltol:11,initial_stol:11,inorgan:[10,11,12],input:2,insert:10,instal:8,instanc:[10,11],instead:[1,2,11],institut:13,integ:10,integr:[1,6,7],intercal:[11,13],intercalant_formula:[4,11],interest:10,intermedi:[1,10,11],internet:11,interpenetr:13,interrupt:8,is_dimension:13,is_intercal:13,is_interpenetr:13,is_vdw_heterostructur:10,isomorph:11,isort:10,issu:2,istructur:11,iter:11,its:[2,6,7,8,11],iupac:[1,10,11,12],jain:[3,6,7,13],json:[1,4],just:11,kei:[1,10,11,12,13],keyword:[10,13],known:[10,11,12],lab:3,label:[1,10,11,12,13],laboratori:8,lanthanid:[10,11,12],largest:11,last:[1,11],latex:[1,10],latexifi:[1,10],latter:11,lawrenc:[3,8],lbl:[10,11,13],lbnl:13,led:3,length:[1,6,7,10,11,12],less:[10,12],levi:11,liabil:8,liabl:8,librari:11,licens:2,life:[6,7],like:[10,11,12,13],likeness_tol:11,likes:11,limit:[8,11],line:[6,7,10],list:[1,2,6,7,8,10,11,12,13],load:10,load_condensed_structure_json:10,local:[6,7,11],local_env:[10,11],logo:[6,7],look:[6,7,12],loss:8,love:2,mai:[8,11],main:[2,10],maintain:[2,6,7],make:8,makeup:[1,10,12],manag:[6,7],mani:[1,2,6,7],map:[10,11],markup:10,mass:11,master:2,match:[1,10,11,12,13],match_bond_angl:11,match_bond_dist:11,match_n_sp:11,materi:[1,3,6,7,8,11],materialsproject:1,matmin:[11,13],max_n_match:11,maximum:11,mean:[2,11],meaning:2,measur:11,mehl:11,merchant:8,met:8,method:[10,11],might:2,miner:[1,4,9,10,12],mineral_dict:12,mineral_match:[10,11],mineral_name_constraint:11,mineralmatch:[10,11],minim:[6,7],minimum:11,minimum_geometry_op:11,miscellan:10,mixtur:[6,7],mmc:4,mnm:[6,7],mo4:[4,11],modif:8,modifi:[6,7,8],modul:[6,7,9],molecul:[1,6,7,9,10,12,13],molecular:11,molecular_compon:11,molecule_graph:11,molecule_graph_to_smil:11,molecule_nam:[4,12,13],moleculegraph:11,moleculenam:11,molecules_first:11,molybdenit:4,moment:2,more:[2,6,7,10,11,12,13],mos2:[4,10,11],most:1,move:2,mpid:1,mprester:[6,7],mrc:[6,7],multipl:[10,11],must:[8,11],my_structur:[6,7],mystructur:1,n_species_type_match:4,n_species_types_match:11,name:[6,7,8,10,11,12,13],name_prefer:11,name_sourc:11,nation:[3,8],ndarrai:11,near_neighbor:[10,11],nearest:[1,10,11,12,13],nearneighbor:[10,11],necessari:2,need:2,neg:11,neglig:8,neighbor:[1,10,11,12,13],neighborsitedetail:12,neighbour:11,neither:8,nest:10,next:[1,10,11,12],nextneighborsitedetail:12,nn_sites_a:11,nn_sites_b:11,nn_summaries_match:11,nnn:[4,11],nnn_distanc:4,nnn_sites_a:11,nnn_sites_b:11,nnn_summaries_match:11,nomenclatur:[10,11,12],non:[8,11],none:[4,6,7,10,11,12,13],nor:8,note:[1,2,11],notic:8,noun:[6,7],nsite:12,num_neighbor:13,num_repetit:[4,11],number:[1,10,11,12,13],numoi:11,numpi:11,obj:10,object:[1,6,7,10,11,12],oblig:8,obtain:[1,10,11,13],occur:[10,11],octahedr:[6,7,10,11,12,13],octahedra:[6,7,10,12],onc:11,one:[1,6,7,10,11,12],onli:[1,10,11,12,13],only_describe_bonds_onc:[10,12],only_describe_cation_polyhedra_connect:[10,12],op_index:11,open:[2,8],openbabel:[6,7],option:[6,7,10,11,12,13],order:[1,10,11,12],order_paramet:11,ordered_compon:11,ore:10,org:[1,6,7],orient:[4,6,7,10,12],other:[2,6,7,8,13],other_compon:11,otherwis:8,our:[2,6,7],out:8,output:[11,13],overal:[10,11],overrid:11,owner:8,oxi:1,oxid:[1,10,11],p4_2:[6,7],p6_3:4,packag:[1,6,7,9],page:[1,6,7],paramet:[1,10,11,12,13],part:[10,11,12],particular:8,pass:[10,11,13],patch:8,pentagon:[4,11],percentag:13,perfect:[10,11,12],perform:[8,10],period:[10,11,12],permiss:8,permit:8,perpetu:8,pip:[6,7],place:[1,6,7],planar:[6,7,11],pleas:[2,6,7],plural:10,poly_formula:[4,11,12],polyhedr:12,polyhedra:[1,10,12,13],poscar:[1,6,7],posit:11,possibl:[2,8,10,11],pre:11,precis:1,precomput:11,preferenti:[10,11],prepar:8,present:10,preset:11,primari:3,primarili:3,prior:8,procedur:2,process:1,procur:8,produc:[1,10,11,12,13],product:8,profit:8,program:3,progress:2,project:[1,6,7],promot:8,properti:[1,10,11,13],propos:2,prototyp:[10,11,12],prototype_match:11,provid:[1,2,6,7,8,11,12],pubchem:11,publicli:8,pull:12,purpos:8,push:2,put:11,pymatgen:[1,6,7,10,11],pyramid:[4,11],python:[1,2],question:2,quick:2,quicker:11,rang:1,raw:[1,10,12],read:[6,7],readi:1,real:[6,7],receipt:8,recommend:[6,7,10,11,12],reconstruct:11,recurs:10,redistribut:8,reduc:[10,11],reduced_formula:10,refer:[10,11,12,13],regent:8,releas:[6,7],relev:[10,12],reli:11,remov:10,repeat:11,repeating_unit:[4,11],repetit:11,repo:2,repres:[10,12],represent:[1,6,7,10,11],reproduc:[2,8],requir:[6,7,8,10,13],research:3,reserv:8,resolv:[10,12,13],resourc:2,respect:[1,11],respons:2,result:[1,6,7],retain:8,return_part:[10,12],review:2,right:8,robocri:[6,7],robocrysfeatur:13,robocrystallograph:[1,3,8,10,11,13],row:[10,11,12],royalti:8,rst:[10,11],run:[1,2,6,7],rutil:[6,7],s828:11,same:[1,2,6,7,10,11,12],sampl:2,scienc:11,script:10,second:[1,11],section:1,see:[1,2,10,11],self:11,separ:[8,11],seri:[10,11,12],servic:[2,8],set:[1,10,11,12,13],shall:8,shape:10,share:[6,7,10,11,12,13],should:[2,10,11,12,13],similar:[6,7,10,11,12],simpli:[1,6,7],simplifi:[1,4,11,12],simplify_molecul:[10,11],singl:[10,11,12],single_compon:[10,12],site:[1,4,9,10,12,13],site_index:[10,11,12],site_indic:[11,12],site_summari:11,siteanalyz:11,sitegroup:12,sitestatsfingerprint:11,six:[6,7],skip_fil:10,smallest:11,smile:11,sn2:10,sno2:[6,7,10],sno6:[6,7],sno:[1,6,7,10],softwar:8,some:[2,11],somewan:[6,7],sort:11,sourc:[2,8,10,11,12,13],space:[2,6,7,10,11,12],spacegroup:10,spacegroup_symbol:10,spacegroupanalyz:11,speci:[10,11],special:[8,10],specif:[2,8,11,13],specifi:[1,6,7,13],spg_analyz:11,spg_symbol:[4,10],split:1,squar:11,stack:[2,11],standard:1,stat:11,state:[1,2,10,11],std_dev:11,step:[1,2],still:[2,6,7],store:11,str:[10,11,12,13],strict:8,string:[10,11,12,13],structur:[1,6,7,10,11,12,13],structure_a:11,structure_b:11,structure_match:11,structurecondens:[6,7,10,11,12,13],structuredescrib:[6,7,10,12],structuregraph:11,structurematch:11,style:[2,8],sub:10,subject:8,sublicens:8,submit:2,submodul:9,subpackag:9,subscript:10,substitut:8,summari:[2,10,11,12],superscript:10,superscript_numb:10,support:[1,6,7,10,11],sym_label:[4,10,11,12],symbol:[10,12,13],symmetr:[10,11],symmetri:[1,6,7,10,11,12],symmetry_label:11,symprec:[1,10,11],system:[6,7,10,11],tab:2,tabl:[10,11,12],take:[11,13],team:2,templat:2,test:2,tetragon:[6,7],text:[6,7,10,12],than:[10,11,12],thei:11,theori:8,thereof:8,thi:[1,2,6,7,8,10,11,12,13],think:2,thorough:2,three:[6,7],through:[1,6,7,8,11],tilt:[6,7,10,12,13],time:[2,11],tip:2,to_sit:[10,11,12],togeth:[11,12],toher:11,tol:1,toler:[1,10,11],tool:[1,6,7,11],tort:8,total:12,trace:2,track:[2,6,7],tradit:11,transform:[10,11],transpar:2,treat:11,tri:[2,11],trigon:[6,7],tupl:[11,12],turn:[6,7,11],two:[1,6,7,10,11,12],type:[1,4,10,11,12,13],uncom:[6,7],under:[2,6,7,8,13],understand:2,unicod:[1,10],unicodeifi:10,unicodeify_spacegroup:10,union:[10,11,12,13],unit:11,unittest:2,univers:8,until:11,updat:[1,2],upgrad:8,use:[1,2,6,7,8,10,11],use_common_formula:[10,11],use_conventional_cel:[10,11],use_fingerprint_match:11,use_iupac_formula:[10,11],use_iupac_ord:12,use_online_pubchem:11,use_oxi_st:10,use_structure_graph:11,use_sym_label:10,use_symmetry_equivalent_sit:[10,11],use_symmetry_inequivalnt_sit:11,use_weight:11,used:[1,6,7,8,10,11],using:[1,3,6,7,10,11,13],util:[6,7,9,11,12],val:11,valu:[10,11,12,13],van:11,varieti:10,vdw:[10,11],vdw_heterostructure_info:4,verbos:1,version:1,voronoinn:[10,11],waal:11,wai:[2,8],want:2,warranti:8,websit:11,weight:11,welcom:[2,6,7],well:2,what:2,whatsoev:8,when:[1,2,6,7,10,11,12],where:[2,6,7,11],whether:[1,2,8,10,11,12,13],which:[1,10,11,13],why:2,within:11,without:[8,11],work:[2,6,7,8,11],would:[1,2,6,7,11],write:2,written:8,you:[2,6,7,8],your:[2,6,7,8],zero:11},titles:["&lt;no title&gt;","<code class=\"docutils literal notranslate\"><span class=\"pre\">robocrys</span></code> program","Contributing to robocrystallographer","Contributors","Condensed structure format","Index","Introduction","Introduction","License","robocrys","robocrys package","robocrys.condense package","robocrys.describe package","robocrys.featurize package"],titleterms:{"new":[6,7],acknowledg:[6,7],adapt:[10,12,13],addit:2,argument:1,basic:1,bug:2,cite:[6,7],cli:10,code:2,command:1,compon:11,condens:[1,4,11],content:[1,10,11,12,13],contribut:[2,6,7],contributor:3,describ:[1,12],discuss:2,exampl:[6,7],featur:13,fingerprint:11,format:[4,6,7],get:2,github:2,great:2,help:2,how:[2,6,7],index:5,instal:[6,7],interfac:[1,6,7],intermedi:[6,7],introduct:[6,7],json:[6,7],licens:[6,7,8],line:1,make:2,miner:11,modif:2,modul:[10,11,12,13],molecul:11,name:1,option:1,output:[6,7],packag:[10,11,12,13],posit:1,program:1,pull:2,python:[6,7],refer:2,report:2,request:2,robocri:[1,9,10,11,12,13],robocrystallograph:[2,6,7],site:11,structur:4,submodul:[10,11,12,13],subpackag:10,tabl:1,through:2,todo:[11,12],usag:[1,6,7],util:10,what:[6,7]}})