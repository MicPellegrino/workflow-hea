from prepare import PrepareDropletSubstrate

ifd="/home/michele/workflow-hea/ti-al-droplets/droplets/020/molten_AlTi_droplet_eq.data"
ifs="/media/michele/LaCie/BackupClusterENI/mdmc-slabs/WMo/eam-zhou/BCC-small/110/WMo_110_eq.data"

PrepareDropletSubstrate(ifd,ifs,out_file_droplet="droplet_test.data",out_file_surface="surface_test.data",gap_z=0)