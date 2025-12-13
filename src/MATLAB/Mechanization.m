
% programmer: Arvin chamasemani
classdef Mechanization
    methods(Static)
        % In Mechanization class (add this method)
        function [r_next, v_next, q_next, yaw_next, pitch_next, roll_next] = step(r, v, q, f_b, w_b_ib, dt)
            [r_next, v_next, q_next] = Mechanization.dynamics(r, v, q, f_b, w_b_ib, dt);
            [yaw_next, pitch_next, roll_next] = CoordTransform.q2euler(q_next);
        end

        function [r_n_cal, v_n_cal, q_cal, yaw_cal, pitch_cal, roll_cal] = ...
                 runBatch(t, r_n, v_n, yaw, pitch, roll, f_b, w_b_ib )
            N=length(t);
            r_n_0=r_n(:,1);
            v_n_0=v_n(:,1);
            q_0=CoordTransform.euler2q(yaw(1),pitch(1),roll(1));
            
            r_n_cal(:,1)=r_n_0;
            v_n_cal(:,1)=v_n_0;
            q_cal(:,1)=q_0;

            
            for i=2:N
                dt_i=t(5)-t(4);
                [r_n_cal(:,i),v_n_cal(:,i),q_cal(:,i)]=Mechanization.dynamics( ...
                    r_n_cal(:,i-1),v_n_cal(:,i-1),q_cal(:,i-1),f_b(:,i-1), ...
                    w_b_ib(:,i-1),dt_i);
            end

            yaw_cal = zeros(1, N);
            pitch_cal=zeros(1, N);
            roll_cal= zeros(1, N);
            for i = 1:N
                [yaw_cal(1,i),pitch_cal(1,i),roll_cal(1,i)] = CoordTransform.q2euler(q_cal(:,i));   
            end
           

        end


        %%%%%%
        function  [r_n,v_n,q]=dynamics(r_n,v_n,q,f_b,omega_ib_b,dt)
            % calculating next step quaternion
            omega_e=7.2921158e-5;
            [M,N]=Mechanization.calM_N(r_n(1));
            omega_ie_n = [omega_e*cos(r_n(1)), 0, -omega_e*sin(r_n(1))]';
            omega_en_n=[v_n(2)/(N+r_n(3)), -v_n(1)/(M+r_n(3)), -v_n(2)*tan(r_n(1))/(N+r_n(3))]';
            dteta_ib_b=omega_ib_b*dt;

            C_nb=CoordTransform.q2mat(q);
            C_bn=C_nb';
            dteta_nb_b = dteta_ib_b - C_bn*(omega_ie_n + omega_en_n)*dt;
            dteta=sqrt(dteta_nb_b(1)^2+dteta_nb_b(2)^2+dteta_nb_b(3)^2);

            s = 2*sin(dteta/2)/dteta; c = 2*(cos(dteta/2)-1);

            temp=[c, s*dteta_nb_b(3), -s*dteta_nb_b(2), s*dteta_nb_b(1)
               -s*dteta_nb_b(3), c, s*dteta_nb_b(1), s*dteta_nb_b(2)
               s*dteta_nb_b(2), -s*dteta_nb_b(1), c, s*dteta_nb_b(3)
               -s*dteta_nb_b(1), -s*dteta_nb_b(2), -s*dteta_nb_b(3), c];

            q=q + 0.5*temp*q;
            q = q / norm(q);

            % calculating next step v_n
            temp=[1, dteta_nb_b(3)/2, -dteta_nb_b(2)/2
                -dteta_nb_b(3)/2, 1, dteta_nb_b(1)/2
                dteta_nb_b(2)/2, -dteta_nb_b(1)/2, 1];
            dv_fb = f_b*dt;
            dv_fn = C_nb*temp*dv_fb;

            a1 = 9.7803267715;a2 = 0.0052790414;a3 = 0.0000232718;
            a4 =-0.0000030876910891;a5 = 0.0000000043977311;a6 = 0.0000000000007211;

            gamma=a1*(1+a2*sin(r_n(1))^2+a3*sin(r_n(1))^4)+ (a4+a5*sin(r_n(1))^2)*r_n(3)+a6*r_n(3)^2;
            gamma_n=[0,0,gamma]';
            dv_n= dv_fn-Mechanization.SkewSymmetry(2*omega_ie_n+ omega_en_n)*v_n*dt+gamma_n*dt;
            v_n_old=v_n;
            
            v_n=v_n+dv_n;

            % calculating next step r_n
            temp=[1/(M+r_n(3)), 0, 0;
               0, 1/(N+r_n(3))/cos(r_n(1)), 0 
               0, 0, -1];
            
            r_n = r_n + 0.5*temp*(v_n_old + v_n)*dt;
        end
        function [M,N] = calM_N( phi )
            %IERS (2003)
            a=6378136.6;
            b=6356751.9;
            e=(1-(b/a)^2)^.5;
            N=a/(1-(e*sin(phi))^2)^0.5;
            M=a*(1-e^2)/(1-(e*sin(phi))^2)^1.5;
        end
        function S=SkewSymmetry(v)
            %--------------------------------------------------------------------------
            % SkewSymmetry  Generates a 3×3 skew-symmetric matrix from a 3×1 vector.
            %
            %   S = SkewSymmetry(v)
            %
            %   This function creates the skew-symmetric matrix corresponding to
            %   a 3-element vector v = [v1; v2; v3]. The skew-symmetric matrix S
            %   is defined such that for any 3×1 vector a:
            %
            %       S * a = cross(v, a)
            %
            %   This operator is widely used in rotational kinematics and dynamics,
            %   for example in quaternion and angular velocity relationships:
            %
            %       ω× = SkewSymmetry(ω)
            %
            %--------------------------------------------------------------------------
            % INPUT:
            %   v - 3×1 vector [v1; v2; v3]
            %
            % OUTPUT:
            %   S - 3×3 skew-symmetric matrix
            %       S = [  0   -v3   v2
            %              v3    0  -v1
            %             -v2   v1    0 ]
            %
            S=[0 , -v(3) , v(2);
                v(3) , 0 , -v(1);
                -v(2) , v(1) , 0];
        end
    end
end

