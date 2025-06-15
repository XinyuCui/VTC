function [e2e_expected_delay_of_packet_of_proposed,packet_loss_of_packet_of_proposed,hops] = e2e_expected_delay_and_packet_loss_of_packet(number_of_node,number_of_video,cache_state,T,route_policy,G_with_propogation_delay,G_with_packet_loss,G_with_relay_process_delay,maximum_times_of_ARQ,FEC_parameter,bandwidth_of_transmitter,packet_size,route_policy_flag)
%根据缓存状态、网络拓扑图图的边属性生成各node在packet level的delay和packet loss开始
%route_policy可以是G_with_packet_loss，也可以是G_with_propogation_delay
    e2e_expected_delay_of_packet_of_proposed=zeros(number_of_node,number_of_video);%用于记录各个卫星传所缓存的各个视频（大概率只缓存了一部分）分别需要多少的端到端时延
    packet_loss_of_packet_of_proposed=zeros(number_of_node,1);
    hops=zeros(number_of_node,1);

    for node_index=1:1:number_of_node
        %发送方S的ID为node_index，接收方D的ID为T+2，生成发送方到接收方的最短属性路径
        send=node_index;
        receive=T+2;

        [minimum_value,tree]=dijkstra(route_policy,send);

        route=[];
        u=receive;
        while(u~=send) 
            route=[u route]; 
            u=tree(u); 
        end
        route=[send route];

        %统计总跳数、各跳的时延、各跳的丢包率
        number_of_hop=size(route,2)-1;
        delay_of_each_hop=zeros(number_of_hop,1);
        packet_loss_of_each_hop=zeros(number_of_hop,1);
        for hop_index=1:1:number_of_hop
            delay_of_each_hop(hop_index)=G_with_propogation_delay(route(hop_index),route(hop_index+1));
            packet_loss_of_each_hop(hop_index)=G_with_packet_loss(route(hop_index),route(hop_index+1));
            relay_process_delay_of_each_hop=G_with_relay_process_delay(route(hop_index),route(hop_index+1));
        end

        number_of_delivered_packet=1;
        %number_of_delivered_packet=cache_state(LEO_index,content_index)/packet_size;    
        
        for video_index=1:1:number_of_video
            if cache_state(node_index,video_index)==0
                e2e_expected_delay_of_packet_of_proposed(node_index,video_index)=nan;
            else
                if route_policy_flag==1
                    [e2e_expected_delay_of_packet_of_proposed(node_index,video_index),packet_loss_of_packet_of_proposed(node_index),hops(node_index)]=e2e_delay_of_specified_route_V2(number_of_hop,maximum_times_of_ARQ,delay_of_each_hop,packet_loss_of_each_hop,relay_process_delay_of_each_hop,FEC_parameter,bandwidth_of_transmitter,number_of_delivered_packet,packet_size);
                else
                    [e2e_expected_delay_of_packet_of_proposed(node_index,video_index),packet_loss_of_packet_of_proposed(node_index),hops(node_index)]=e2e_delay_of_specified_route_V2_distance_first(number_of_hop,maximum_times_of_ARQ,delay_of_each_hop,packet_loss_of_each_hop,relay_process_delay_of_each_hop,FEC_parameter,bandwidth_of_transmitter,number_of_delivered_packet,packet_size);
                end
            end
        end
    end        
end

