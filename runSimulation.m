function [ray, meas, MEAS0] = runSimulation(SystemParam, ray, meas, Bounds, Global_Index, Tally, ...
                                    gg_max, I_enter, Theta_enter, ...
                                    aa, xx, yy, h, y0)
    for gg=1:gg_max
        %initialize the entering light values
        I_ent=I_enter(gg);
        Theta_ent=Theta_enter(gg);
        %preallocate empty vectors with limits for the  while loop
        max_loop=1+(SystemParam.maxBounce*SystemParam.scatterNum);%max number of loops to run
        
        I_in=zeros(1,(max_loop));
        I_intrack=zeros(1,(max_loop));
        P_in=zeros((max_loop),2);
        Theta_in=zeros(1,(max_loop));
        Direction_in=zeros(1,(max_loop));
        %assign starting while loop values
        ray.Pow_enter(aa,xx,yy)=I_ent+ray.Pow_enter(aa,xx,yy);
        I_in(1,1)=I_ent;
        Theta_in(1,1)=Theta_ent;
        P_in(1,1:2)=[0,y0(xx,yy,h)];
        Direction_in(1,1)=1;
        %assign initial index values
        while_num=1;
        bounce_num=1;
        st_index=2;%starting index to allocate the first scattered ray we're tracking
        pow_check=0;
        MEAS0=sum(meas.inten);
        while any(I_in) && while_num<=max_loop
            %%%%%%5%initial variable set up for while loop%%%%%
            I_0=I_in(1,while_num);
            I_intrack(1,while_num)=I_0;
            pow_check=pow_check+I_0;
            P_0=P_in(while_num,:);
            Theta_0=Theta_in(1,while_num);
            direction=Direction_in(1,while_num);
            V_0=[direction*abs(cos(Theta_0)),sin(Theta_0)];
            %initial measurement for use in later
            %troubleshooting/difference checking
            
            meas_start0=sum(meas.inten);
            b2h0=ray.b2hpow(aa,xx,yy);
            cutoff0=ray.cutoffpow(aa,xx,yy);
            trans0=ray.transmitted(aa,xx,yy);
            abs0=ray.absorbed(aa,xx,yy);
            back0=ray.backscat(aa,xx,yy);
            smaabs0=ray.SMAabs(aa,xx,yy);
            approx0=ray.approxpow(aa,xx,yy);
            if I_0<SystemParam.intensityMin%check to make sure the while loop is worth it to keep running
                ray.cutoffpow(aa,xx,yy)=ray.cutoffpow(aa,xx,yy)+I_0;%mark amount of light left after cutoff
                I_in(1,while_num)=0;%redundant with the continue but just in case
                while_num=while_num+1;%update the while number
                continue%go to the next while loop
            elseif length(I_0)>1 || length((P_0))>2 || length(Theta_0)>1 || length(V_0)>2
                
                error('too big vector')
            elseif isnan(I_0)
                
                error('is NaN')
                
            end
            
            %%%%%%traveling %%%%%%%%%%%%
            %traveling the length of the fiber
            %start by recording the initial currently measured
            %intensity
            meas_start1=sum(meas.inten);
            %                     %check transmission %should be in the travel function
            [I_1,V_1,P_1,direction,meas,IT,Tally] = Traveling(I_0,V_0,P_0,direction,SystemParam,Bounds,Global_Index,Tally,meas);%temporary
            %record various loss tallies
            ray.SMAabs(aa,xx,yy)=IT.housi+ray.SMAabs(aa,xx,yy);
            ray.absorbed(aa,xx,yy)=IT.absorbi+ray.absorbed(aa,xx,yy);
            ray.backscat(aa,xx,yy)=IT.backi+ray.backscat(aa,xx,yy);
            ray.cutoffpow(aa,xx,yy)=IT.cutoffi+ray.cutoffpow(aa,xx,yy);
            ray.approxpow(aa,xx,yy)=IT.approxi+ray.approxpow(aa,xx,yy);
            ray.b2hpow(aa,xx,yy)=IT.b2hi+ray.b2hpow(aa,xx,yy);
            ray.transmitted(aa,xx,yy)=IT.transi+ray.transmitted(aa,xx,yy);
            %after travel function
            meas_dif=sum(meas.inten)-meas_start1;
            %                     disp('Traveling')
            Pow_out=[meas_dif,I_1,IT.backi,IT.absorbi,IT.housi,IT.cutoffi,IT.approxi];
            
            % [Tally,Differenceamount,Diffamountpos] =  DifTrack(I_0,Pow_out,SystemParam,'Traveling',0,Global_Index,Tally);
            % ray.approxpow_dif(aa,xx,yy)=Differenceamount+ray.approxpow_dif(aa,xx,yy);%sum of power differences due to difference in approximations
            % ray.approxpow_pos(aa,xx,yy)=Diffamountpos+ray.approxpow_pos(aa,xx,yy);%sum of total power diff
            
            %check if there's enough power to warrant going thru the end reflect
            %function
            if I_1<SystemParam.intensityMin
                %store and record data. much more probably
                ray.cutoffpow(aa,xx,yy)=ray.cutoffpow(aa,xx,yy)+I_1;%
                %lots of stuff in here
                %reset loop
                I_in(1,while_num)=0;
                while_num=while_num+1;
                
                continue%go to the next while loop
            end
            meas_dif_1=sum(meas.inten)-meas_start0;
            b2h_1=ray.b2hpow(aa,xx,yy)-b2h0;
            cutoff_1=ray.cutoffpow(aa,xx,yy)-cutoff0;
            trans_1=ray.transmitted(aa,xx,yy)-trans0;
            abs_1=ray.absorbed(aa,xx,yy)-abs0;
            smaabs_1=ray.SMAabs(aa,xx,yy)-smaabs0;
            approx_1=ray.approxpow(aa,xx,yy)-approx0;
            back_1=ray.backscat(aa,xx,yy)-back0;
            %                     disp('while loop step 0 to step 1')
            Pow_Out=[meas_dif_1,I_1,b2h_1,cutoff_1,trans_1,abs_1,smaabs_1,approx_1,back_1];
            %                     approx_dif_final=ray.approxpow_dif(aa,xx,yy)-approxdif0;
    % %                                 [Tally,~,~] =  DifTrack(I_0,Pow_Out,SystemParam,'while loop step 0 to step 1',0,Global_Index,Tally);
            
            cutoffstart=ray.cutoffpow(aa,xx,yy);
            meas_start2=sum(meas.inten);%initial total measured before entering into the end reflect function
            %%%%%%%end reflect stuff%%%%%%%%%
            %               disp('sum forward scatter')
            if bounce_num<=SystemParam.maxBounce
                %create the new rays for the loop, update the
                %bounce number, take relevant measurements and
                %tallies.
                [I_scat,Theta_scat,I_trans,direction,b2h,bounce_num,Tally,meas]= EndReflect(I_1,V_1,P_1,SystemParam,direction,bounce_num,Global_Index,Tally,meas);
                %record and save
                ray.transmitted(aa,xx,yy)=ray.transmitted(aa,xx,yy)+I_trans;
                ray.b2hpow(aa,xx,yy)=b2h+ray.b2hpow(aa,xx,yy);%record and save
                I_scat(I_scat<SystemParam.intensityMin)=0;%any of the scattered rays less than the minimum tracking power are set to 0
                %create approprate indices for assigning the scattered values to the main while loop vectors
                end_index=st_index+length(I_scat)-1;
                indices_assign=st_index:1:end_index;
                indices_assign=indices_assign(indices_assign<=max_loop);%make sure  the indices are under the max value
                %sort the incoming vector to prioritize highest
                %bvalues to be assigned to the vector
                [I_sort_scat,sort_scat]=sort(I_scat,'descend');
                
                if ~isempty(indices_assign)%otherwise, nothing will need to be assigned to the loop vector
                    if end_index>max_loop%if there aren't enough spots to allocate all new I_sorted
                        highest_index=find(indices_assign==max_loop);%gives the total amount of spots available to allocate the sorted scI_Scat to
                        ray.cutoffpow(aa,xx,yy)=sum(I_scat(sort_scat(highest_index+1:end)))+ray.cutoffpow(aa,xx,yy);
                        sort_scat=sort_scat(1:highest_index);
                    else%if theres plenty of space
                        highest_index=length(indices_assign);
                    end
                    
                    %new loop assigned rays
                    I_in(1,indices_assign(1):indices_assign(highest_index))=I_scat(sort_scat);
                    Theta_in(1,indices_assign(1):indices_assign(highest_index))=Theta_scat(sort_scat);
                    P_in(indices_assign(1):indices_assign(highest_index),1)=P_1(1);P_in(indices_assign(1):indices_assign(end),2)=P_1(2);%all will have same initial point
                    Direction_in(1,indices_assign(1):indices_assign(highest_index))=direction;
                    sum_forwardscat=sum(I_scat(sort_scat));
                else
                    sum_forwardscat=0;
                    
                end
                %update the starting index
                st_index=end_index+1;
            else%if we've already had the max number of bounces, just record the transmitted and/or lost light
                [I_scat,~,I_trans,~,b2h,~,Tally,meas]= EndReflect(I_1,V_1,P_1,SystemParam,direction,bounce_num,Global_Index,Tally,meas);
                ray.transmitted(aa,xx,yy)=ray.transmitted(aa,xx,yy)+I_trans;
                ray.b2hpow(aa,xx,yy)=b2h+ray.b2hpow(aa,xx,yy);%record and save
                %record what isn't being tracked
                ray.cutoffpow(aa,xx,yy)=sum(I_scat)+ray.cutoffpow(aa,xx,yy);
                sum_forwardscat=0;
            end
            %check if there's any major differences in the total
            meas_dif=sum(meas.inten)-meas_start2;
            cut_off=ray.cutoffpow(aa,xx,yy)-cutoffstart;
            %                     disp('end reflect')
            Pow_Out=[meas_dif,sum_forwardscat,cut_off,I_trans,b2h];
            dif_endref=I_1-sum(Pow_Out);
    % %                                 [Tally,Differenceamount,Diffamountpos] =  DifTrack(I_1,Pow_Out,SystemParam,'End Reflect',0,Global_Index,Tally);
    % %                                 ray.approxpow_dif(aa,xx,yy)=Differenceamount+ray.approxpow_dif(aa,xx,yy);%sum of power differences due to difference in approximations
    % %                                 ray.approxpow_pos(aa,xx,yy)=Diffamountpos+ray.approxpow_pos(aa,xx,yy);%sum of total power diff
            
            %set the just used I_in to 0
            I_in(1,while_num)=0;
            
            %update the while counter
            while_num=while_num+1;
            
            if while_num<(max_loop)%more than one while loop remaining
                %sort all of the remaining in descending order so the
                %highest intensity ones are prioritized
                [I_in(1,while_num:end),sort_Index]=sort(I_in(while_num:end),'descend');
                sort_Index=sort_Index+(while_num-1);
                Theta_in(1,while_num:end)=Theta_in(1,sort_Index);
                P_in(while_num:end,1)=P_in(sort_Index,1);P_in(while_num:end,2)=P_in(sort_Index,2);
                Direction_in(1,while_num:end)=Direction_in(1,sort_Index);
            end
            
            %if the loop cuts off with anything remaining in it
            if while_num==(2+(SystemParam.maxBounce*SystemParam.scatterNum)) && any(I_in)%these conditions shouldn't occur BUT
                ray.cutoffpow(aa,xx,yy)=sum(I_in)+ray.cutoffpow(aa,xx,yy);
            end
            
            meas_dif_final=sum(meas.inten)-meas_start0;
            b2h_final=ray.b2hpow(aa,xx,yy)-b2h0;
            cutoff_final=ray.cutoffpow(aa,xx,yy)-cutoff0;
            trans_final=ray.transmitted(aa,xx,yy)-trans0;
            abs_final=ray.absorbed(aa,xx,yy)-abs0;
            smaabs_final=ray.SMAabs(aa,xx,yy)-smaabs0;
            approx_final=ray.approxpow(aa,xx,yy)-approx0;
            back_final=ray.backscat(aa,xx,yy)-back0;
            %                     disp('while loop')
            %                     start=[I_0,meas_start0,b2h0,cutoff0,trans0,abs0,smaabs0,approx0]
            Pow_out=[sum_forwardscat,meas_dif_final,b2h_final,cutoff_final,trans_final,abs_final,smaabs_final,approx_final,back_final];
            %                     approx_dif_final=ray.approxpow_dif(aa,xx,yy)-approxdif0;
    % %                                 [Tally,~,~] =  DifTrack(I_0,Pow_Out,SystemParam,'while loop',0,Global_Index,Tally);
            
        end
        %summing all of the measured side emitted power
    end
end
