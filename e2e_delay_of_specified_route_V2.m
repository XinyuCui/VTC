function [expected_e2e_delay,packet_loss_with_HARQ,hops] = e2e_delay_of_specified_route_V2(number_of_hop,maximum_times_of_ARQ,delay_of_each_hop,packet_loss_of_each_hop,relay_process_delay_of_each_hop,FEC_parameter,bandwidth_of_transmitter,number_of_delivered_packet,packet_size)
%计算未使用HARQ前，接收方的累积丢包率
packet_not_loss_without_HARQ=1;
for i=1:1:number_of_hop
    packet_not_loss_without_HARQ=packet_not_loss_without_HARQ*(1-packet_loss_of_each_hop(i));
end
packet_loss_without_HARQ=1-packet_not_loss_without_HARQ;

%计算使用FEC后，接收方的丢包率
packet_loss_with_FEC=(1-FEC_parameter)*packet_loss_without_HARQ;

%计算在FEC的基础上使用ARQ后，即HARQ后，接收方的丢包率
packet_not_loss_with_HARQ=0;
for i=0:1:maximum_times_of_ARQ
    packet_not_loss_with_HARQ=packet_not_loss_with_HARQ+(1-packet_loss_with_FEC)*packet_loss_with_FEC^i;
end
packet_loss_with_HARQ=1-packet_not_loss_with_HARQ;

%计算在HARQ机制下，packet端到端的时延
accumulated_one_way_delay=sum(delay_of_each_hop);
accumulated_relay_process_delay=sum(relay_process_delay_of_each_hop);
expected_transmission_times_of_packet=1;

for i=1:1:maximum_times_of_ARQ
    expected_transmission_times_of_packet=expected_transmission_times_of_packet+i*(1-packet_loss_with_FEC)*packet_loss_with_FEC^(i);
end

%expected_e2e_delay=packet_size*number_of_delivered_packet*expected_transmission_times_of_packet/bandwidth_of_transmitter+accumulated_one_way_delay+accumulated_relay_process_delay;%这么做相当于每次等用户接收到包后再发下一个包，是不对的
expected_e2e_delay=packet_size*number_of_delivered_packet*expected_transmission_times_of_packet/bandwidth_of_transmitter;%端到端时延的期望值

hops=number_of_hop;
end

