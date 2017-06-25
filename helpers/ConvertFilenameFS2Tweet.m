function res_str = ConvertFilenameFS2Tweet(str_orig,bird_name,serial_num)
Y = str_orig(1:4);
M = str_orig(6:7);
D = str_orig(9:10);
h = str_orig(12:13);
m = str_orig(15:16);
s = str_orig(18:19);
res_str = [bird_name '_' num2str(serial_num) '_' Y '_' M '_' D '_' h '_' m '_' s];