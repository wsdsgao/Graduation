% function [de_out] = GMSK_demod_new(signal_recv, de_out, n, viterbi_deep, g)
    %myFun - Description
    %
    % Syntax: output = myFun(input)
    %
    % Long description
    clc;
    clear all; 
    close all;

    signal_recv = [0.923880911075620 - 0.382680104199425i,0.707131559550881 - 0.707082001953902i,0.382980125856272 - 0.923756582222348i,0.00225129112370721 - 0.999997465840927i,-0.371624142171022 - 0.928383270506127i,-0.672059635718298 - 0.740497026353373i,-0.855014345944814 - 0.518604346519158i,-0.928055858252627 - 0.372441034209417i,-0.925978043405724 - 0.377577360458101i,-0.847164384406339 - 0.531330881648554i,-0.656733916081332 - 0.754122379636408i,-0.351861476621691 - 0.936052082562505i,0.0146171781915893 - 0.999893163343822i,0.360793496621042 - 0.932645727377744i,0.611107442078465 - 0.791547657590063i,0.731835717947539 - 0.681481094335133i,0.728046892724362 - 0.685527331325602i,0.599154778923734 - 0.800633218704327i,0.341030512594971 - 0.940052226995506i,-0.0108375059064624 - 0.999941272508405i,-0.384534322484795 - 0.923110694787554i,-0.707302019773008 - 0.706911488678055i,-0.923888993447641 - 0.382660590845601i,-0.999999999997107 - 2.40522719106862e-06i,-0.923880862418421 + 0.382680221669469i,-0.707109280452249 + 0.707104281912013i,-0.382686698382020 + 0.923878179676016i,-3.60784078561082e-06 + 0.999999999993492i,0.382679110597507 + 0.923881322634190i,0.707103431532267 + 0.707110130824960i,0.923877719451122 + 0.382687809450203i,0.999999999988430 + 4.81045438095269e-06i,0.923881782848623 - 0.382677999524992i,0.707110981196649 - 0.707102581151500i,0.382688920517832 - 0.923877259224891i,6.01306797576518e-06 - 0.999999999981922i,-0.382676888451923 - 0.923882243061720i,-0.707101730769709 - 0.707111831567314i,-0.923876798997325 - 0.382690031584906i,-0.999999999973967 - 7.21568157073060e-06i];
    de_out = zeros(1,3040);
    n = 11;
    viterbi_deep = 5;
    load('lib/g_1024.mat');  % GMSK调制 g函数 
    g = g(1:16:end);

    init_path = zeros(8, 3);
    path_record = cell(32, 1);
    path_record_tmp = cell(32, 1);
    path_weight = zeros(32, viterbi_deep);
    
    init_path(1,:) = [0,0,0];
    init_path(2,:) = [0,0,1];
    init_path(3,:) = [0,1,0];
    init_path(4,:) = [0,1,1];
    init_path(5,:) = [1,0,0];
    init_path(6,:) = [1,0,1];
    init_path(7,:) = [1,1,0];
    init_path(8,:) = [1,1,1];

    if n == 1
        init_path = [repmat(zeros(1,2),8,1), init_path];
    else
        init_path = [repmat(de_out(n-2:n-1),8,1), init_path];
    end

    % 初始化放入path_record中

    phi_last = 0;
    for j = 1:8
        path_record{bin2dec(num2str(init_path(j,:)))+1} = init_path(j,:);
    end

    path_record_tmp = path_record;

    for i = 1:viterbi_deep
        signal_recv_ang = angle(signal_recv(1+(i-1)*4:i*4));
        for j = 1 : 32
            if ~isempty(path_record{j})
                bit_record = path_record{j};
                bit_5 = 2*bit_record - 1;
                g_path_ang_dif = zeros(1,2);
                signal_recv_ang_dif = zeros(1,2);
                if i == 1
                    [~, I_sig, Q_sig, ~] = GMSK(bit_5, phi_last, g);
                    g_path_ang = angle(complex(I_sig, Q_sig));
                    % figure;
                    % plot(ang)
                    % hold on;
                    % plot(signal_recv_ang)
                    % hold off;
                
                    for m = 1:2
                        g_path_ang_dif(m) = g_path_ang(m) - g_path_ang(m+2);
                        if  g_path_ang_dif(m) > pi
                            g_path_ang_dif(m) =  g_path_ang_dif(m) - 2*pi;
                        elseif  g_path_ang_dif(m) < -pi
                            g_path_dif_dif(m) =  g_path_ang_dif(m) + 2*pi;
                        end
                        signal_recv_ang_dif(m) = signal_recv_ang(m) - signal_recv_ang(m+2);
                        if signal_recv_ang_dif(m) > pi
                            signal_recv_ang_dif(m) = signal_recv_ang_dif(m) - 2*pi;
                        elseif signal_recv_ang_dif(m) < -pi
                            signal_recv_ang_dif(m) = signal_recv_ang_dif(m) + 2*pi;
                        end
                    end

                    path_weight(j,i) =  path_weight(j,i) + (sum(signal_recv_ang_dif) - sum(g_path_ang_dif)).^2;
                else
                    bit_record_cur0 = [bit_record, 0];
                    bit_record_cur1 = [bit_record, 1];
                    [~, I_sig0, Q_sig0, ~] = GMSK([bit_5(end-3:end),-1], 0, g);
                    [~, I_sig1, Q_sig1, ~] = GMSK([bit_5(end-3:end),1], 0, g);
                    g_path_ang0 = angle(complex(I_sig0, Q_sig0));
                    g_path_ang1 = angle(complex(I_sig1, Q_sig1));
                    state_index0 = bin2dec(num2str([bit_record(end-3:end),0]))+1;
                    state_index1 = bin2dec(num2str([bit_record(end-3:end),1]))+1;

                    for m = 1:2
                        g_path_ang_dif0(m) = g_path_ang0(m) - g_path_ang0(m+2);
                        if  g_path_ang_dif0(m) > pi
                            g_path_ang_dif0(m) =  g_path_ang_dif0(m) - 2*pi;
                        elseif  g_path_ang_dif0(m) < -pi
                            g_path_dif_dif0(m) =  g_path_ang_dif0(m) + 2*pi;
                        end
                        g_path_ang_dif1(m) = g_path_ang1(m) - g_path_ang1(m+2);
                        if  g_path_ang_dif1(m) > pi
                            g_path_ang_dif1(m) =  g_path_ang_dif1(m) - 2*pi;
                        elseif  g_path_ang_dif1(m) < -pi
                            g_path_dif_dif1(m) =  g_path_ang_dif1(m) + 2*pi;
                        end
                        signal_recv_ang_dif(m) = signal_recv_ang(m) - signal_recv_ang(m+2);
                        if signal_recv_ang_dif(m) > pi
                            signal_recv_ang_dif(m) = signal_recv_ang_dif(m) - 2*pi;
                        elseif signal_recv_ang_dif(m) < -pi
                            signal_recv_ang_dif(m) = signal_recv_ang_dif(m) + 2*pi;
                        end
                    end

                    if path_weight(state_index0,i) == 0 || path_weight(j,i-1) + (sum(signal_recv_ang_dif) - sum(g_path_ang_dif0)).^2 < path_weight(state_index0,i)
                        path_weight(state_index0,i) =  path_weight(j,i-1) + (sum(signal_recv_ang_dif) - sum(g_path_ang_dif0)).^2;
                        path_record_tmp{state_index0} = bit_record_cur0;                   
                    end
                    
                    if path_weight(state_index1,i) == 0 || path_weight(j,i-1) + (sum(signal_recv_ang_dif) - sum(g_path_ang_dif1)).^2 < path_weight(state_index1,i)
                        path_weight(state_index1,i) =  path_weight(j,i-1) + (sum(signal_recv_ang_dif) - sum(g_path_ang_dif1)).^2;
                        path_record_tmp{state_index1} = bit_record_cur1;                   
                    end
                end
            end
        end

        path_record = path_record_tmp;
    end

    [~,index_min] = min(path_weight(:,end))
    out = path_record{index_min}
    de_out(n:n+4) = out(3:end-2);
    path_record
% end