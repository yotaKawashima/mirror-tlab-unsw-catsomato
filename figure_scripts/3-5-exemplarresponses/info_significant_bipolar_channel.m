% Try all possible bipolar channels
s1_f1_23_logsnr = [32, 53, 59, 60, 61, 89, 122, ...
                   131, 133, 137, 141, 154, 173, 174, 175, 176];
               
s2_f1_23_logsnr = [2, 7, 60, 63, 85];

s1_f1_23_velogp = [23, 32, 39, 49, 53, 60, 61, 69, 71, 121,...
                   139, 141, 142, 154, 163, 164, 169, 175, 176, 177];
               
s2_f1_23_velogp = [2, 11, 15, 17, 60, 63, 67, 74, 76, 85];

s1_f2_200_logsnr = [26, 27, 32, 36, 45, 50, 54, 60, 61, 63, 69, 72, 77, ...
                    121, 122, 123, 129, 130, 132, 137, 141, 143, 149, ...
                    150, 151, 155, 159, 160, 165, 168, 177, ...
                    178, 180];
s2_f2_200_logsnr = [4, 12 19, 21, 25, 31, 36, 42, 43, 44, 47, 48, 49,...
                    54, 60, 61, 62, 63, 69, 70, 74, 77, 82, 83, 86,...
                    87, 91, 92, 93, 95, 97, 99, 100, 103, 104, 106, 112];

s1_f2_200_velogp = [23, 24, 25, 32, 39, 40, 41, 45, 48, 52, 54, 60, ...
                    61, 63, 69, 71, 72, 77, 86, 87, 121, 122, 129, ...
                    130, 132, 137, 141, 143, 149, 150, 151, 159, ...
                    160, 168, 171, 177, 178];

s2_f2_200_velogp = [2, 4, 12, 19, 21, 23, 24, 25,31, 36, 44, 46, 47, ...
                    48, 49, 61, 62, 63, 69, 70, 74, 77, 82, 83, 84, ...
                    86, 87, 91, 92, 93, 99, 100, 104];


plotpair = {};
increment = 1;
freqs = [23, 200];

for i_type = 1:2
    for i_freq = 1:2
        for i_area = 1:2
            switch i_type*4 + i_freq*2 + i_area
                case 7 % (1, 1, 1)
                    channel_set = s1_f1_23_logsnr;
                case 8 % (1, 1, 2)
                    channel_set = s2_f1_23_logsnr;
                case 9 % (1, 2, 1)
                    channel_set = s1_f2_200_logsnr;
                case 10 % (1, 2, 2)
                    channel_set = s2_f2_200_logsnr;
                case 11 % (2, 1, 1)
                    channel_set = s1_f1_23_velogp;
                case 12 % (2, 1, 2)
                    channel_set = s2_f1_23_velogp;
                case 13 % (2, 2, 1)
                    channel_set = s1_f2_200_velogp;
                case 14 % (2, 2, 2)
                    channel_set = s2_f2_200_velogp;                    
            end
            
            for i_channel = 1:length(channel_set)
                bipolar_channel = channel_set(i_channel);
                freq = freqs(i_freq);
                plotpair{increment} = [i_area, bipolar_channel, ...
                                       freq, i_type];
                increment = increment + 1; 
            end % i_channel = 1:length(
        end % i_freq = 1:2
    end % i_area =1:2
end % i_type = 1:2

clear increment 
clear channel_set
clear freq
clear bipolar_channel
clear freqs
