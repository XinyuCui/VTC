%V3是在V2的基础上，将1000个用户分别请求1~1000个视频，改为了100个用户，随机的请求视频
clc;clear all;clf;
close all;

%adjacency_matrix_of_each_time_step(1,:,:)是一个（1，68，68）维的数组
%在用的时候注意用squeeze函数对维度压缩一下，压缩成（68,68）维的，squeeze函数的作用就是去除数组中维度为1的部分
%同理sat_positions也是
load('adjacency_matrix_of_each_time_step.mat')
load('sat_positions_of_each_time_step.mat')

%系统参数设置
%T、total_time、time_step取决于satellite_topology_graph_data_generating.m这个程序里是怎么设的
T = 66; % 卫星总数
total_time = 240; % 总模拟时间 (s)
time_step = 2; % 时间步长 (s)
num_steps = total_time / time_step;

maximum_times_of_ARQ=4;%最大重传次数
FEC_parameter=0.2;%FEC可以恢复x%的丢包
bandwidth_of_transmitter=10;%单位Mbps
packet_size=0.001;%单位Mb
number_of_video=1000;%视频总数
size_of_video=1000;%视频大小1000Mb
cache_capacity=number_of_video*size_of_video*0.1;%单个活动基站的缓存容量为：可缓存远端视频库中10%的视频
number_of_node=T+1;%即T个LEO+1个远端视频库
number_of_LEO=T;
number_of_request=100;%这里假设每个用户只发起一个请求

%*****************************Begin:用户请求生成阶段*************************************
%暂时不用基于流行度的请求生成，采用随机的请求生成
video_request=randi([1,number_of_video],number_of_request,1);
%*****************************Stop:用户请求生成阶段*************************************

%*****************************Begin:缓存阶段*************************************
%由于分发阶段缓存状态是不变的，所以在蒙特卡洛仿真前就把缓存放置好
%**********暂时采用随机缓存策略进行缓存放置**********
cache_state=zeros(number_of_node,number_of_video);
cache_state(end,:)=size_of_video;%备注，最后一行是远端视频库，存了所有的视频，其他行才是低轨卫星存的

for node_index=1:1:number_of_node-1
   random_video_ID=randperm(number_of_video,number_of_video);%生成指定范围内的随机数
   for video_index=1:1:number_of_video
       cache_state(node_index,random_video_ID(video_index))=rand*size_of_video;

       if sum(cache_state(node_index,:))>cache_capacity
           cache_state(node_index,random_video_ID(video_index))=cache_state(node_index,random_video_ID(video_index))-(sum(cache_state(node_index,:))-cache_capacity);%这一步是为了保证分配给各视频的缓存容量之和正好等于cache capacity
           break;
       end
   end
end

for node_index=1:1:number_of_node
   for video_index=1:1:number_of_video
       if cache_state(node_index,video_index)<1
           cache_state(node_index,video_index)=0;%按上述方法弄出来的cache_state,有些是特别小的负数
       end
       
       if isnan(cache_state(node_index,video_index))
           cache_state(node_index,video_index)=0;
       end
   end
end
%*****************************Stop:缓存阶段*************************************

%蒙特卡洛仿真设置
monte_carlo_iteration=60;

monte_carlo_average_delay_of_proposed=zeros(monte_carlo_iteration,1);
monte_carlo_average_delay_of_comparison=zeros(monte_carlo_iteration,1);
monte_carlo_average_delay_of_random_dis=zeros(monte_carlo_iteration,1);%random_dis指的是随机选择delivery nodes，relay nodes以distance优先的原则用dijkstra获得
monte_carlo_average_delay_of_random_pac=zeros(monte_carlo_iteration,1);%random_dis指的是随机选择delivery nodes，relay nodes以packet loss优先的原则用dijkstra获得

monte_carlo_average_packet_loss_of_proposed=zeros(monte_carlo_iteration,1);
monte_carlo_average_packet_loss_of_comparison=zeros(monte_carlo_iteration,1);
monte_carlo_average_packet_loss_of_random_dis=zeros(monte_carlo_iteration,1);
monte_carlo_average_packet_loss_of_random_pac=zeros(monte_carlo_iteration,1);

monte_carlo_average_packet_loss_of_proposed_of_hop_level=zeros(number_of_node,monte_carlo_iteration);
monte_carlo_average_packet_loss_of_comparison_of_hop_level=zeros(number_of_node,monte_carlo_iteration);
monte_carlo_average_packet_loss_of_random_dis_of_hop_level=zeros(number_of_node,monte_carlo_iteration);
monte_carlo_average_packet_loss_of_random_pac_of_hop_level=zeros(number_of_node,monte_carlo_iteration);


for monte_carlo_index=1:1:monte_carlo_iteration
    fprintf('%d\n',monte_carlo_index)
    
    %重置缓存状态  
    cache_state_of_proposed=zeros(number_of_request,number_of_node,number_of_video);
    cache_state_of_comparison=zeros(number_of_request,number_of_node,number_of_video);
    cache_state_of_random_pac=zeros(number_of_request,number_of_node,number_of_video);
    cache_state_of_random_dis=zeros(number_of_request,number_of_node,number_of_video); 
    
    for request_index=1:1:number_of_request
        cache_state_of_proposed(request_index,:,:)=cache_state;
        cache_state_of_comparison(request_index,:,:)=cache_state;
        cache_state_of_random_pac(request_index,:,:)=cache_state;
        cache_state_of_random_dis(request_index,:,:)=cache_state;      
    end

    %**********创建视频各项指标的统计数组**********
    e2e_expected_delay_of_each_request_of_proposed=zeros(number_of_request,1);
    e2e_expected_delay_of_each_request_of_comparison=zeros(number_of_request,1);
    e2e_expected_delay_of_each_request_of_random_dis=zeros(number_of_request,1);
    e2e_expected_delay_of_each_request_of_random_pac=zeros(number_of_request,1);
    
    packet_loss_of_each_request_of_proposed=zeros(number_of_request,3);%第一列是丢的包总量，第2列是发的包总量，第3列用于最后统计丢包率
    packet_loss_of_each_request_of_comparison=zeros(number_of_request,3);
    packet_loss_of_each_request_of_random_dis=zeros(number_of_request,3);
    packet_loss_of_each_request_of_random_pac=zeros(number_of_request,3);
    
    packet_loss_of_each_request_of_proposed_of_hop_level=zeros(number_of_node,3);%hop_level的丢包率统计，第一列是丢的包总量，第2列是发的包总量，第3列用于最后统计丢包率
    packet_loss_of_each_request_of_comparison_of_hop_level=zeros(number_of_node,3);
    packet_loss_of_each_request_of_random_dis_of_hop_level=zeros(number_of_node,3);
    packet_loss_of_each_request_of_random_pac_of_hop_level=zeros(number_of_node,3);
    
    user_received_size_of_proposed=zeros(number_of_request,1);        
    user_received_size_of_comparison=zeros(number_of_request,1);     
    user_received_size_of_random_dis=zeros(number_of_request,1);     
    user_received_size_of_random_pac=zeros(number_of_request,1); 
    
    % 循环计算每个时间步
    for step = 1:num_steps
        %首先给每个视频分配时长为time_step的传输时间
        left_time_of_proposed=time_step*ones(number_of_request,1);
        left_time_of_comparison=time_step*ones(number_of_request,1);
        left_time_of_random_dis=time_step*ones(number_of_request,1);
        left_time_of_random_pac=time_step*ones(number_of_request,1);
        
        %调取网络拓扑
        adjacency_matrix=squeeze(adjacency_matrix_of_each_time_step(step,:,:));
        sat_positions=squeeze(sat_positions_of_each_time_step(step,:,:));
        
        
        %*****************************Begin:根据网络拓扑图，生成图的边属性*************************************
        G_with_propogation_delay=zeros(size(adjacency_matrix,1),size(adjacency_matrix,2));
        G_with_packet_loss=zeros(size(adjacency_matrix,1),size(adjacency_matrix,2));
        G_with_relay_process_delay=zeros(size(adjacency_matrix,1),size(adjacency_matrix,2));

        for row_index=1:1:size(adjacency_matrix,1)
           for column_index=1:1:size(adjacency_matrix,2)
               if adjacency_matrix(row_index,column_index)==1
                   G_with_propogation_delay(row_index,column_index)=norm(sat_positions(row_index,:)-sat_positions(column_index,:))/300000;
                   G_with_packet_loss(row_index,column_index)=0.5*rand;%设ISL的packet loss服从0~0.5之间的均匀分布
                   G_with_relay_process_delay(row_index,column_index)=0.001;%设中继处理时延为0.001s
               end
           end
        end
        %*****************************Stop:根据网络拓扑图，生成图的边属性*************************************

        
        %*****************************Begin:根据缓存状态、网络拓扑图图的边属性生成各node在packet level的delay和packet loss*************************************
        %路由按照G_with_packet_loss生成的，去计算各个低轨卫星传各个视频的packet需要多长时间、各个低轨卫星的丢包率
        [e2e_expected_delay_of_packet_of_proposed,packet_loss_of_packet_of_proposed,number_of_hop] = e2e_expected_delay_and_packet_loss_of_packet(number_of_node,number_of_video,cache_state,T,G_with_packet_loss,G_with_propogation_delay,G_with_packet_loss,G_with_relay_process_delay,maximum_times_of_ARQ,FEC_parameter,bandwidth_of_transmitter,packet_size,1);
        %*****************************Stop:根据缓存状态、网络拓扑图图的边属性生成各node在packet level的delay和packet loss*************************************
        
        %*****************************Begin：测试所提方法的性能*****************************
        %**********计算传各个视频需要多长时间、传各个视频会丢百分之多少的包**********
        for request_index=1:1:number_of_request
            [sorted_result,sorted_node_index]=sort(e2e_expected_delay_of_packet_of_proposed(:,video_request(request_index)),'ascend');
            for node_index=1:1:number_of_node
                %首先判断node是否缓存了video，没有的话直接换下一个node
                if isnan(e2e_expected_delay_of_packet_of_proposed(sorted_node_index(node_index),video_request(request_index)))
                    continue;
                end
                
                %首先判断缓存的内容是否足以恢复视频
                if user_received_size_of_proposed(request_index)+cache_state_of_proposed(request_index,sorted_node_index(node_index),video_request(request_index))<size_of_video
                    deliver_time=min(left_time_of_proposed(request_index),cache_state_of_proposed(request_index,sorted_node_index(node_index),video_request(request_index))/packet_size*e2e_expected_delay_of_packet_of_proposed(sorted_node_index(node_index),video_request(request_index)));
                else
                    deliver_time=min(left_time_of_proposed(request_index),(size_of_video-user_received_size_of_proposed(request_index))/packet_size*e2e_expected_delay_of_packet_of_proposed(sorted_node_index(node_index),video_request(request_index)));
                end
                
                %计算传了多少包，丢了多少包。min函数是考虑到剩余时间截止时，有些包可能传不完
                number_of_delivered_packet=deliver_time/e2e_expected_delay_of_packet_of_proposed(sorted_node_index(node_index),video_request(request_index));
                number_of_effective_delivered_packet=number_of_delivered_packet*(1-packet_loss_of_packet_of_proposed(sorted_node_index(node_index)));
                number_of_loss_packet=number_of_delivered_packet*packet_loss_of_packet_of_proposed(sorted_node_index(node_index));
                
                %对缓存中未发的视频大小、用户收到的视频大小、剩余时间进行更新
                cache_state_of_proposed(request_index,sorted_node_index(node_index),video_request(request_index))=cache_state_of_proposed(request_index,sorted_node_index(node_index),video_request(request_index))-number_of_effective_delivered_packet*packet_size;                
                user_received_size_of_proposed(request_index)=user_received_size_of_proposed(request_index)+number_of_effective_delivered_packet*packet_size;
                left_time_of_proposed(request_index)=left_time_of_proposed(request_index)-deliver_time;                
                
                %进行KPI的记录
                e2e_expected_delay_of_each_request_of_proposed(request_index)=e2e_expected_delay_of_each_request_of_proposed(request_index)+deliver_time;                
                packet_loss_of_each_request_of_proposed(request_index,1)=packet_loss_of_each_request_of_proposed(request_index,1)+number_of_loss_packet;
                packet_loss_of_each_request_of_proposed(request_index,2)=packet_loss_of_each_request_of_proposed(request_index,2)+number_of_delivered_packet;
                packet_loss_of_each_request_of_proposed_of_hop_level(number_of_hop(sorted_node_index(node_index)),1)=packet_loss_of_each_request_of_proposed_of_hop_level(number_of_hop(sorted_node_index(node_index)),1)+number_of_loss_packet;
                packet_loss_of_each_request_of_proposed_of_hop_level(number_of_hop(sorted_node_index(node_index)),2)=packet_loss_of_each_request_of_proposed_of_hop_level(number_of_hop(sorted_node_index(node_index)),2)+number_of_delivered_packet;
                
                %如果视频video_index的传输时长time_step耗尽，则结束对视频video_index的传输
                if left_time_of_proposed(request_index)<0.001
                    left_time_of_proposed(request_index)=0;
                    break;
                end
                
                %如果用户已经收到了完整视频video_index，则结束对视频video_index的传输
                if user_received_size_of_proposed(request_index)>=size_of_video
                    break;
                end
            end
        end
        %*****************************Stop：测试所提方法的性能*****************************
           
        

        
               
        
        
        %*****************************Begin:根据缓存状态、网络拓扑图图的边属性生成各node在packet level的delay和packet loss*************************************
        %路由按照G_with_propogation_delay生成的，去计算各个低轨卫星传各个视频的packet需要多长时间、各个低轨卫星的丢包率
        %e2e_expected_delay_of_packet_of_comparison_fake和real分别指的是不考虑丢包的时延和考虑丢包的时延，e2e_expected_delay_and_packet_loss_of_packet函数的最后一个变量分别赋的2和1
        [e2e_expected_delay_of_packet_of_comparison_fake,packet_loss_of_packet_of_comparison,number_of_hop] = e2e_expected_delay_and_packet_loss_of_packet(number_of_node,number_of_video,cache_state,T,G_with_propogation_delay,G_with_propogation_delay,G_with_packet_loss,G_with_relay_process_delay,maximum_times_of_ARQ,FEC_parameter,bandwidth_of_transmitter,packet_size,2);
        [e2e_expected_delay_of_packet_of_comparison_real,packet_loss_of_packet_of_comparison,number_of_hop] = e2e_expected_delay_and_packet_loss_of_packet(number_of_node,number_of_video,cache_state,T,G_with_propogation_delay,G_with_propogation_delay,G_with_packet_loss,G_with_relay_process_delay,maximum_times_of_ARQ,FEC_parameter,bandwidth_of_transmitter,packet_size,1);
        %*****************************Stop:根据缓存状态、网络拓扑图图的边属性生成各node在packet level的delay和packet loss*************************************

        %*****************************Begin：测试对比方法的性能*****************************
        %**********计算传各个视频需要多长时间、传各个视频会丢百分之多少的包**********
        for request_index=1:1:number_of_request
            [sorted_result,sorted_node_index]=sort(e2e_expected_delay_of_packet_of_comparison_fake(:,video_request(request_index)),'ascend');%注意这里用的是e2e_expected_delay_of_packet_of_comparison_fake
            for node_index=1:1:number_of_node
                %首先判断node是否缓存了video，没有的话直接换下一个node
                if isnan(e2e_expected_delay_of_packet_of_comparison_real(sorted_node_index(node_index),video_request(request_index)))
                    continue;
                end
                
                %首先判断缓存的内容是否足以恢复视频
                if user_received_size_of_comparison(request_index)+cache_state_of_comparison(request_index,sorted_node_index(node_index),video_request(request_index))<size_of_video
                    deliver_time=min(left_time_of_comparison(request_index),cache_state_of_comparison(request_index,sorted_node_index(node_index),video_request(request_index))/packet_size*e2e_expected_delay_of_packet_of_comparison_real(sorted_node_index(node_index),video_request(request_index)));
                else
                    deliver_time=min(left_time_of_comparison(request_index),(size_of_video-user_received_size_of_comparison(request_index))/packet_size*e2e_expected_delay_of_packet_of_comparison_real(sorted_node_index(node_index),video_request(request_index)));
                end
                
                %计算传了多少包，丢了多少包。min函数是考虑到剩余时间截止时，有些包可能传不完
                number_of_delivered_packet=deliver_time/e2e_expected_delay_of_packet_of_comparison_real(sorted_node_index(node_index),video_request(request_index));
                number_of_effective_delivered_packet=number_of_delivered_packet*(1-packet_loss_of_packet_of_comparison(sorted_node_index(node_index)));
                number_of_loss_packet=number_of_delivered_packet*packet_loss_of_packet_of_comparison(sorted_node_index(node_index));
                
                %对缓存中未发的视频大小、用户收到的视频大小、剩余时间进行更新
                cache_state_of_comparison(request_index,sorted_node_index(node_index),video_request(request_index))=cache_state_of_comparison(request_index,sorted_node_index(node_index),video_request(request_index))-number_of_effective_delivered_packet*packet_size;                
                user_received_size_of_comparison(request_index)=user_received_size_of_comparison(request_index)+number_of_effective_delivered_packet*packet_size;
                left_time_of_comparison(request_index)=left_time_of_comparison(request_index)-deliver_time;                
                
                %进行KPI的记录
                e2e_expected_delay_of_each_request_of_comparison(request_index)=e2e_expected_delay_of_each_request_of_comparison(request_index)+deliver_time;                
                packet_loss_of_each_request_of_comparison(request_index,1)=packet_loss_of_each_request_of_comparison(request_index,1)+number_of_loss_packet;
                packet_loss_of_each_request_of_comparison(request_index,2)=packet_loss_of_each_request_of_comparison(request_index,2)+number_of_delivered_packet;
                packet_loss_of_each_request_of_comparison_of_hop_level(number_of_hop(sorted_node_index(node_index)),1)=packet_loss_of_each_request_of_comparison_of_hop_level(number_of_hop(sorted_node_index(node_index)),1)+number_of_loss_packet;
                packet_loss_of_each_request_of_comparison_of_hop_level(number_of_hop(sorted_node_index(node_index)),2)=packet_loss_of_each_request_of_comparison_of_hop_level(number_of_hop(sorted_node_index(node_index)),2)+number_of_delivered_packet;

                %如果视频video_index的传输时长time_step耗尽，则结束对视频video_index的传输
                if left_time_of_comparison(request_index)<0.001
                    left_time_of_comparison(request_index)=0;
                    break;
                end
                
                %如果用户已经收到了完整视频video_index，则结束对视频video_index的传输
                if user_received_size_of_comparison(request_index)>=size_of_video
                    break;
                end
            end
        end
        %*****************************Stop：测试对比方法的性能*****************************        


        
        
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
                
        %*****************************Begin:根据缓存状态、网络拓扑图图的边属性生成各node在packet level的delay和packet loss*************************************
        %路由按照G_with_packet_loss生成的，去计算各个低轨卫星传各个视频的packet需要多长时间、各个低轨卫星的丢包率
        [e2e_expected_delay_of_packet_of_random_pac,packet_loss_of_packet_of_random_pac,number_of_hop] = e2e_expected_delay_and_packet_loss_of_packet(number_of_node,number_of_video,cache_state,T,G_with_packet_loss,G_with_propogation_delay,G_with_packet_loss,G_with_relay_process_delay,maximum_times_of_ARQ,FEC_parameter,bandwidth_of_transmitter,packet_size,1);
        %*****************************Stop:根据缓存状态、网络拓扑图图的边属性生成各node在packet level的delay和packet loss*************************************
        
        %*****************************Begin：测试random选择delivery node方法的性能，relay nodes按照packet loss优先的原则选*****************************
        %**********计算传各个视频需要多长时间、传各个视频会丢百分之多少的包**********
        for request_index=1:1:number_of_request
            sorted_node_index=randperm(number_of_node);
            for node_index=1:1:number_of_node
                %首先判断node是否缓存了video，没有的话直接换下一个node
                if isnan(e2e_expected_delay_of_packet_of_random_pac(sorted_node_index(node_index),video_request(request_index)))
                    continue;
                end
                                
                %首先判断缓存的内容是否足以恢复视频
                if user_received_size_of_random_pac(request_index)+cache_state_of_random_pac(request_index,sorted_node_index(node_index),video_request(request_index))<size_of_video
                    deliver_time=min(left_time_of_random_pac(request_index),cache_state_of_random_pac(request_index,sorted_node_index(node_index),video_request(request_index))/packet_size*e2e_expected_delay_of_packet_of_random_pac(sorted_node_index(node_index),video_request(request_index)));
                else
                    deliver_time=min(left_time_of_random_pac(request_index),(size_of_video-user_received_size_of_random_pac(request_index))/packet_size*e2e_expected_delay_of_packet_of_random_pac(sorted_node_index(node_index),video_request(request_index)));
                end
                
                %计算传了多少包，丢了多少包。min函数是考虑到剩余时间截止时，有些包可能传不完
                number_of_delivered_packet=deliver_time/e2e_expected_delay_of_packet_of_random_pac(sorted_node_index(node_index),video_request(request_index));
                number_of_effective_delivered_packet=number_of_delivered_packet*(1-packet_loss_of_packet_of_random_pac(sorted_node_index(node_index)));
                number_of_loss_packet=number_of_delivered_packet*packet_loss_of_packet_of_random_pac(sorted_node_index(node_index));
                
                %对缓存中未发的视频大小、用户收到的视频大小、剩余时间进行更新
                cache_state_of_random_pac(request_index,sorted_node_index(node_index),video_request(request_index))=cache_state_of_random_pac(request_index,sorted_node_index(node_index),video_request(request_index))-number_of_effective_delivered_packet*packet_size;                
                user_received_size_of_random_pac(request_index)=user_received_size_of_random_pac(request_index)+number_of_effective_delivered_packet*packet_size;
                left_time_of_random_pac(request_index)=left_time_of_random_pac(request_index)-deliver_time;  
               
                %进行KPI的记录
                e2e_expected_delay_of_each_request_of_random_pac(request_index)=e2e_expected_delay_of_each_request_of_random_pac(request_index)+deliver_time;                
                packet_loss_of_each_request_of_random_pac(request_index,1)=packet_loss_of_each_request_of_random_pac(request_index,1)+number_of_loss_packet;
                packet_loss_of_each_request_of_random_pac(request_index,2)=packet_loss_of_each_request_of_random_pac(request_index,2)+number_of_delivered_packet;
                packet_loss_of_each_request_of_random_pac_of_hop_level(number_of_hop(sorted_node_index(node_index)),1)=packet_loss_of_each_request_of_random_pac_of_hop_level(number_of_hop(sorted_node_index(node_index)),1)+number_of_loss_packet;
                packet_loss_of_each_request_of_random_pac_of_hop_level(number_of_hop(sorted_node_index(node_index)),2)=packet_loss_of_each_request_of_random_pac_of_hop_level(number_of_hop(sorted_node_index(node_index)),2)+number_of_delivered_packet;

                %如果视频video_index的传输时长time_step耗尽，则结束对视频video_index的传输
                if left_time_of_random_pac(request_index)<0.001
                    left_time_of_random_pac(request_index)=0;
                    break;
                end
                
                %如果用户已经收到了完整视频video_index，则结束对视频video_index的传输
                if user_received_size_of_random_pac(request_index)>=size_of_video
                    break;
                end
            end
        end
        %*****************************Stop：测试random选择delivery node方法的性能，relay nodes按照packet loss优先的原则选*****************************        

        
        
        
        
        
        
        
        
        %*****************************Begin:根据缓存状态、网络拓扑图图的边属性生成各node在packet level的delay和packet loss*************************************
        %路由按照G_with_propogation_delay生成的，去计算各个低轨卫星传各个视频的packet需要多长时间、各个低轨卫星的丢包率
        %e2e_expected_delay_of_packet_of_comparison_fake和real分别指的是不考虑丢包的时延和考虑丢包的时延，e2e_expected_delay_and_packet_loss_of_packet函数的最后一个变量分别赋的2和1
        [e2e_expected_delay_of_packet_of_random_dis,packet_loss_of_packet_of_random_dis,number_of_hop] = e2e_expected_delay_and_packet_loss_of_packet(number_of_node,number_of_video,cache_state,T,G_with_propogation_delay,G_with_propogation_delay,G_with_packet_loss,G_with_relay_process_delay,maximum_times_of_ARQ,FEC_parameter,bandwidth_of_transmitter,packet_size,1);
        %*****************************Stop:根据缓存状态、网络拓扑图图的边属性生成各node在packet level的delay和packet loss*************************************

        %*****************************Begin：测试random选择delivery node方法的性能，relay nodes按照distance优先的原则选*****************************
        for request_index=1:1:number_of_request
            sorted_node_index=randperm(number_of_node);
            for node_index=1:1:number_of_node
                %首先判断node是否缓存了video，没有的话直接换下一个node
                if isnan(e2e_expected_delay_of_packet_of_random_dis(sorted_node_index(node_index),video_request(request_index)))
                    continue;
                end
                
                %首先判断缓存的内容是否足以恢复视频
                if user_received_size_of_random_dis(request_index)+cache_state_of_random_dis(request_index,sorted_node_index(node_index),video_request(request_index))<size_of_video
                    deliver_time=min(left_time_of_random_dis(request_index),cache_state_of_random_dis(request_index,sorted_node_index(node_index),video_request(request_index))/packet_size*e2e_expected_delay_of_packet_of_random_dis(sorted_node_index(node_index),video_request(request_index)));
                else
                    deliver_time=min(left_time_of_random_dis(request_index),(size_of_video-user_received_size_of_random_dis(request_index))/packet_size*e2e_expected_delay_of_packet_of_random_dis(sorted_node_index(node_index),video_request(request_index)));
                end                
                
                %计算传了多少包，丢了多少包。min函数是考虑到剩余时间截止时，有些包可能传不完
                number_of_delivered_packet=deliver_time/e2e_expected_delay_of_packet_of_random_dis(sorted_node_index(node_index),video_request(request_index));
                number_of_effective_delivered_packet=number_of_delivered_packet*(1-packet_loss_of_packet_of_random_dis(sorted_node_index(node_index)));
                number_of_loss_packet=number_of_delivered_packet*packet_loss_of_packet_of_random_dis(sorted_node_index(node_index));
                
                %对缓存中未发的视频大小、用户收到的视频大小、剩余时间进行更新
                cache_state_of_random_dis(request_index,sorted_node_index(node_index),video_request(request_index))=cache_state_of_random_dis(request_index,sorted_node_index(node_index),video_request(request_index))-number_of_effective_delivered_packet*packet_size;                
                user_received_size_of_random_dis(request_index)=user_received_size_of_random_dis(request_index)+number_of_effective_delivered_packet*packet_size;
                left_time_of_random_dis(request_index)=left_time_of_random_dis(request_index)-deliver_time;                
                
                %进行KPI的记录
                e2e_expected_delay_of_each_request_of_random_dis(request_index)=e2e_expected_delay_of_each_request_of_random_dis(request_index)+deliver_time;                
                packet_loss_of_each_request_of_random_dis(request_index,1)=packet_loss_of_each_request_of_random_dis(request_index,1)+number_of_loss_packet;
                packet_loss_of_each_request_of_random_dis(request_index,2)=packet_loss_of_each_request_of_random_dis(request_index,2)+number_of_delivered_packet;
                packet_loss_of_each_request_of_random_dis_of_hop_level(number_of_hop(sorted_node_index(node_index)),1)=packet_loss_of_each_request_of_random_dis_of_hop_level(number_of_hop(sorted_node_index(node_index)),1)+number_of_loss_packet;
                packet_loss_of_each_request_of_random_dis_of_hop_level(number_of_hop(sorted_node_index(node_index)),2)=packet_loss_of_each_request_of_random_dis_of_hop_level(number_of_hop(sorted_node_index(node_index)),2)+number_of_delivered_packet;

                %如果视频video_index的传输时长time_step耗尽，则结束对视频video_index的传输
                if left_time_of_random_dis(request_index)<0.001
                    left_time_of_random_dis(request_index)=0;
                    break;
                end
                
                %如果用户已经收到了完整视频video_index，则结束对视频video_index的传输
                if user_received_size_of_random_dis(request_index)>=size_of_video
                    break;
                end
            end
        end
        %*****************************Stop：测试random选择delivery node方法的性能，relay nodes按照distance优先的原则选*****************************
    
    end
    
    %*****************************Begin:时延结果统计*************************************
    average_delay_of_proposed=sum(e2e_expected_delay_of_each_request_of_proposed)/number_of_request;
    average_delay_of_comparison=sum(e2e_expected_delay_of_each_request_of_comparison)/number_of_request;
    average_delay_of_random_dis=sum(e2e_expected_delay_of_each_request_of_random_dis)/number_of_request;
    average_delay_of_random_pac=sum(e2e_expected_delay_of_each_request_of_random_pac)/number_of_request;
    %*****************************Stop:时延结果统计*************************************

    %*****************************Begin:丢包率结果统计*************************************
    for request_index=1:1:number_of_request
        packet_loss_of_each_request_of_proposed(request_index,3)=packet_loss_of_each_request_of_proposed(request_index,1)/packet_loss_of_each_request_of_proposed(request_index,2);
        packet_loss_of_each_request_of_comparison(request_index,3)=packet_loss_of_each_request_of_comparison(request_index,1)/packet_loss_of_each_request_of_comparison(request_index,2);
        packet_loss_of_each_request_of_random_dis(request_index,3)=packet_loss_of_each_request_of_random_dis(request_index,1)/packet_loss_of_each_request_of_random_dis(request_index,2);
        packet_loss_of_each_request_of_random_pac(request_index,3)=packet_loss_of_each_request_of_random_pac(request_index,1)/packet_loss_of_each_request_of_random_pac(request_index,2);
    end
    
    average_packet_loss_of_request_of_proposed=sum(packet_loss_of_each_request_of_proposed(:,3))/number_of_request;
    average_packet_loss_of_request_of_comparison=sum(packet_loss_of_each_request_of_comparison(:,3))/number_of_request;
    average_packet_loss_of_request_of_random_dis=sum(packet_loss_of_each_request_of_random_dis(:,3))/number_of_request;
    average_packet_loss_of_request_of_random_pac=sum(packet_loss_of_each_request_of_random_pac(:,3))/number_of_request;
    
    for hop_index=1:1:number_of_node
        packet_loss_of_each_request_of_proposed_of_hop_level(hop_index,3)=packet_loss_of_each_request_of_proposed_of_hop_level(hop_index,1)/packet_loss_of_each_request_of_proposed_of_hop_level(hop_index,2);
        packet_loss_of_each_request_of_comparison_of_hop_level(hop_index,3)=packet_loss_of_each_request_of_comparison_of_hop_level(hop_index,1)/packet_loss_of_each_request_of_comparison_of_hop_level(hop_index,2);
        packet_loss_of_each_request_of_random_dis_of_hop_level(hop_index,3)=packet_loss_of_each_request_of_random_dis_of_hop_level(hop_index,1)/packet_loss_of_each_request_of_random_dis_of_hop_level(hop_index,2);
        packet_loss_of_each_request_of_random_pac_of_hop_level(hop_index,3)=packet_loss_of_each_request_of_random_pac_of_hop_level(hop_index,1)/packet_loss_of_each_request_of_random_pac_of_hop_level(hop_index,2);
    end   
    %*****************************Stop:丢包率结果统计*************************************

    %*****************************Begin:将时延、丢包率结果加入蒙特卡洛数组*************************************
    monte_carlo_average_delay_of_proposed(monte_carlo_index)=average_delay_of_proposed;
    monte_carlo_average_delay_of_comparison(monte_carlo_index)=average_delay_of_comparison;
    monte_carlo_average_delay_of_random_dis(monte_carlo_index)=average_delay_of_random_dis;
    monte_carlo_average_delay_of_random_pac(monte_carlo_index)=average_delay_of_random_pac;

    monte_carlo_average_packet_loss_of_proposed(monte_carlo_index)=average_packet_loss_of_request_of_proposed;
    monte_carlo_average_packet_loss_of_comparison(monte_carlo_index)=average_packet_loss_of_request_of_comparison;
    monte_carlo_average_packet_loss_of_random_dis(monte_carlo_index)=average_packet_loss_of_request_of_random_dis;
    monte_carlo_average_packet_loss_of_random_pac(monte_carlo_index)=average_packet_loss_of_request_of_random_pac;
    
    for hop_index=1:1:number_of_node
        monte_carlo_average_packet_loss_of_proposed_of_hop_level(hop_index,monte_carlo_index)=packet_loss_of_each_request_of_proposed_of_hop_level(hop_index,3);
        monte_carlo_average_packet_loss_of_comparison_of_hop_level(hop_index,monte_carlo_index)=packet_loss_of_each_request_of_comparison_of_hop_level(hop_index,3);
        monte_carlo_average_packet_loss_of_random_dis_of_hop_level(hop_index,monte_carlo_index)=packet_loss_of_each_request_of_random_dis_of_hop_level(hop_index,3);
        monte_carlo_average_packet_loss_of_random_pac_of_hop_level(hop_index,monte_carlo_index)=packet_loss_of_each_request_of_random_pac_of_hop_level(hop_index,3);
    end   
    %*****************************Stop:将时延、丢包率结果加入蒙特卡洛数组*************************************
end


%*****************************Begin:计算时延、视频丢包率、各跳丢包率的蒙特卡洛均值*************************************
average_monte_carlo_average_delay_of_proposed=sum(monte_carlo_average_delay_of_proposed)/monte_carlo_iteration;
average_monte_carlo_average_delay_of_comparison=sum(monte_carlo_average_delay_of_comparison)/monte_carlo_iteration;
average_monte_carlo_average_delay_of_random_dis=sum(monte_carlo_average_delay_of_random_dis)/monte_carlo_iteration;
average_monte_carlo_average_delay_of_random_pac=sum(monte_carlo_average_delay_of_random_pac)/monte_carlo_iteration;

average_monte_carlo_average_packet_loss_of_proposed=sum(monte_carlo_average_packet_loss_of_proposed)/monte_carlo_iteration;
average_monte_carlo_average_packet_loss_of_comparison=sum(monte_carlo_average_packet_loss_of_comparison)/monte_carlo_iteration;
average_monte_carlo_average_packet_loss_of_random_dis=sum(monte_carlo_average_packet_loss_of_random_dis)/monte_carlo_iteration;
average_monte_carlo_average_packet_loss_of_random_pac=sum(monte_carlo_average_packet_loss_of_random_pac)/monte_carlo_iteration;


average_monte_carlo_average_packet_loss_of_proposed_of_hop_level=zeros(number_of_node,1);
average_monte_carlo_average_packet_loss_of_comparison_of_hop_level=zeros(number_of_node,1);
average_monte_carlo_average_packet_loss_of_random_dis_of_hop_level=zeros(number_of_node,1);
average_monte_carlo_average_packet_loss_of_random_pac_of_hop_level=zeros(number_of_node,1);

for hop_index=1:1:number_of_node
    average_monte_carlo_average_packet_loss_of_proposed_of_hop_level(hop_index)=nanmean(monte_carlo_average_packet_loss_of_proposed_of_hop_level(hop_index,:));
    average_monte_carlo_average_packet_loss_of_comparison_of_hop_level(hop_index)=nanmean(monte_carlo_average_packet_loss_of_comparison_of_hop_level(hop_index,:));
    average_monte_carlo_average_packet_loss_of_random_dis_of_hop_level(hop_index)=nanmean(monte_carlo_average_packet_loss_of_random_dis_of_hop_level(hop_index,:));
    average_monte_carlo_average_packet_loss_of_random_pac_of_hop_level(hop_index)=nanmean(monte_carlo_average_packet_loss_of_random_pac_of_hop_level(hop_index,:));
end
%*****************************Stop:计算时延、视频丢包率、各跳丢包率的蒙特卡洛均值结束*************************************


