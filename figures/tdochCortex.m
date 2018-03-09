function s = tdochCortex(r, pars)

    if nargin == 1, pars = loadParameters(); end
  
    s = initialiseStates(r, pars);
    s = corticalProcessing(r, s, pars);

end



function s = corticalProcessing(r, s, pars)

   if pars.verb, fprintf('Simulating cortical processing...');  end;

    nOfIts = length(r.timeSpace) - 1;
    markerThresh = nOfIts / 20;
    marker = 0;

    for ti = 1:(length(r.timeSpace) - 1)

        if pars.verb
            marker = marker + 1;
            if marker > markerThresh
                fprintf('.');
                marker = 0;
            end
        end

        dt = r.timeSpace(ti + 1) - r.timeSpace(ti);
        s = updateStateVars(s, r, dt, pars);

    end

    if pars.verb, fprintf('..done! <3 time = %.1fm\n', toc / 60); end;

end



function s = initialiseStates(r, pars)

    % Parameters
    outOfClockIterations = 2000;

    dt = r.timeSpace(3) - r.timeSpace(2);

    s0.t = 1;
    s0.n.SPn = zeros(outOfClockIterations + 1, pars.N);
    s0.n.SPa = zeros(outOfClockIterations + 1, pars.N);

    s0.p.Sg  = zeros(outOfClockIterations + 1, pars.N);
    s0.p.Sn  = zeros(outOfClockIterations + 1, pars.N);
    s0.p.Sa  = zeros(outOfClockIterations + 1, pars.N);
    s0.p.Hi  = zeros(outOfClockIterations + 1, pars.N);
    s0.p.He  = zeros(outOfClockIterations + 1, pars.N);
    s0.p.xi  = zeros(outOfClockIterations + 1, pars.N);
    s0.p.xe  = zeros(outOfClockIterations + 1, pars.N);
    s0.p.Ai  = zeros(outOfClockIterations + 1, pars.N);
    s0.p.Ae  = zeros(outOfClockIterations + 1, pars.N);
    
    s0.q.Sg  = zeros(outOfClockIterations + 1, pars.N);
    s0.q.Sn  = zeros(outOfClockIterations + 1, pars.N);
    s0.q.Sa  = zeros(outOfClockIterations + 1, pars.N);
    s0.q.Hi  = zeros(outOfClockIterations + 1, pars.N);
    s0.q.He  = zeros(outOfClockIterations + 1, pars.N);
    s0.q.xi  = zeros(outOfClockIterations + 1, pars.N);
    s0.q.xe  = zeros(outOfClockIterations + 1, pars.N);
    s0.q.Ai  = zeros(outOfClockIterations + 1, pars.N);
    s0.q.Ae  = zeros(outOfClockIterations + 1, pars.N);

    r0 = r;
    r0.A = zeros(outOfClockIterations + 1, pars.N);
    r0.E = zeros(outOfClockIterations + 1, 1);

    for ii = 1:outOfClockIterations
        s0 = updateStateVars(s0, r0, dt, pars);
    end

    s.n.SPn = zeros(length(r.timeSpace), pars.N);
    s.n.SPa = zeros(length(r.timeSpace), pars.N);

    s.p.Sg  = zeros(length(r.timeSpace), pars.N);
    s.p.Sn  = zeros(length(r.timeSpace), pars.N);
    s.p.Sa  = zeros(length(r.timeSpace), pars.N);
    s.p.Hi  = zeros(length(r.timeSpace), pars.N);
    s.p.He  = zeros(length(r.timeSpace), pars.N);
    s.p.xi  = zeros(length(r.timeSpace), pars.N);
    s.p.xe  = zeros(length(r.timeSpace), pars.N);
    s.p.Ai  = zeros(length(r.timeSpace), pars.N);
    s.p.Ae  = zeros(length(r.timeSpace), pars.N);

    s.q.Sg  = zeros(length(r.timeSpace), pars.N);
    s.q.Sn  = zeros(length(r.timeSpace), pars.N);
    s.q.Sa  = zeros(length(r.timeSpace), pars.N);
    s.q.Hi  = zeros(length(r.timeSpace), pars.N);
    s.q.He  = zeros(length(r.timeSpace), pars.N);
    s.q.xi  = zeros(length(r.timeSpace), pars.N);
    s.q.xe  = zeros(length(r.timeSpace), pars.N);
    s.q.Ai  = zeros(length(r.timeSpace), pars.N);
    s.q.Ae  = zeros(length(r.timeSpace), pars.N);

    s.t = 1;
    s.n.SPn(1, :) = s0.n.SPn(end, :); 
    s.n.SPa(1, :) = s0.n.SPa(end, :); 

    s.p.Sg(1, :)  = s0.p.Sg(end, :); 
    s.p.Sn(1, :)  = s0.p.Sn(end, :); 
    s.p.Sa(1, :)  = s0.p.Sa(end, :); 
    s.p.Hi(1, :)  = s0.p.Hi(end, :);
    s.p.He(1, :)  = s0.p.He(end, :);
    s.p.xi(1, :)  = s0.p.xi(end, :);
    s.p.xe(1, :)  = s0.p.xe(end, :);
    s.p.Ai(1, :)  = s0.p.Ai(end, :);
    s.p.Ae(1, :)  = s0.p.Ae(end, :);

    s.q.Sg(1, :)  = s0.q.Sg(end, :); 
    s.q.Sn(1, :)  = s0.q.Sn(end, :); 
    s.q.Sa(1, :)  = s0.q.Sa(end, :); 
    s.q.Hi(1, :)  = s0.q.Hi(end, :);
    s.q.He(1, :)  = s0.q.He(end, :);
    s.q.xi(1, :)  = s0.q.xi(end, :);
    s.q.xe(1, :)  = s0.q.xe(end, :);
    s.q.Ai(1, :)  = s0.q.Ai(end, :);
    s.q.Ae(1, :)  = s0.q.Ae(end, :);

end



function s = updateStateVars(s, r, dt, pars)

	tt = s.t + 1;

    [xPi, xPe, xQi, xQe] = inputCurrent(s, dt, pars);
    s.p.xi(tt, :)  = xPi;
    s.p.xe(tt, :)  = xPe;
    s.q.xi(tt, :)  = xQi;
    s.q.xe(tt, :)  = xQe;

    [dSn, dSp, dSq] = dSdt(s, r, pars);
    s.n.SPn(tt, :) = s.n.SPn(s.t, :) + dSn.Pn * dt;
    s.n.SPa(tt, :) = s.n.SPa(s.t, :) + dSn.Pa * dt;
    s.p.Sg(tt, :)  = s.p.Sg(s.t, :)  + dSp.g  * dt;
    s.p.Sn(tt, :)  = s.p.Sn(s.t, :)  + dSp.n  * dt;
    s.p.Sa(tt, :)  = s.p.Sa(s.t, :)  + dSp.a  * dt;
    s.q.Sg(tt, :)  = s.q.Sg(s.t, :)  + dSq.g  * dt;
    s.q.Sn(tt, :)  = s.q.Sn(s.t, :)  + dSq.n  * dt;
    s.q.Sa(tt, :)  = s.q.Sa(s.t, :)  + dSq.a  * dt;

    % We need to ensure S > 0 because noise can be negative
	s.n.SPn(tt, :) = heaviside(s.n.SPn(tt, :)) .* s.n.SPn(tt, :);
	s.n.SPa(tt, :) = heaviside(s.n.SPa(tt, :)) .* s.n.SPa(tt, :);
	s.p.Sg(tt, :)  = heaviside(s.p.Sg(tt, :))  .* s.p.Sg(tt, :);
	s.p.Sn(tt, :)  = heaviside(s.p.Sn(tt, :))  .* s.p.Sn(tt, :);
	s.p.Sa(tt, :)  = heaviside(s.p.Sa(tt, :))  .* s.p.Sa(tt, :);
	s.q.Sg(tt, :)  = heaviside(s.q.Sg(tt, :))  .* s.q.Sg(tt, :);
	s.q.Sn(tt, :)  = heaviside(s.q.Sn(tt, :))  .* s.q.Sn(tt, :);
	s.q.Sa(tt, :)  = heaviside(s.q.Sa(tt, :))  .* s.q.Sa(tt, :);

    [dAPi, dAPe, dAQi, dAQe] = dAdt(s, pars);
    s.p.Ai(tt, :) = s.p.Ai(s.t, :) + dAPi * dt;
    s.p.Ae(tt, :) = s.p.Ae(s.t, :) + dAPe * dt;
    s.q.Ai(tt, :) = s.q.Ai(s.t, :) + dAQi * dt;
    s.q.Ae(tt, :) = s.q.Ae(s.t, :) + dAQe * dt;

    [dHPi, dHPe, dHQi, dHQe] = dHdt(s, pars);
    s.p.Hi(tt, :) = s.p.Hi(s.t, :) + dHPi * dt;
    s.p.He(tt, :) = s.p.He(s.t, :) + dHPe * dt; 
    s.q.Hi(tt, :) = s.q.Hi(s.t, :) + dHQi * dt;
    s.q.He(tt, :) = s.q.He(s.t, :) + dHQe * dt; 

    s.t = tt;

end



function [xPi, xPe, xQi, xQe] = inputCurrent(s, dt, pars)

    J = 1:pars.N;
    I = 1:pars.N;

    % Delays
    minDelay = pars.selfDelay / dt;
    maxDelay = pars.maxDelay  / dt;
    delay = round(linspace(minDelay, maxDelay, pars.N));
    tij = s.t - delay(abs(repmat(I',size(J)) - repmat(J',size(I))') + 1);
    tij = max(1, tij);
    indPS  = sub2ind(size(s.p.Sn), tij, repmat(J, [pars.N, 1]));
    indQS  = max(1, round(s.t - minDelay));
    
    % Pitch inputs
    % Thalamic-to-pitch
    thalP = pars.Jtn * s.n.SPn(s.t, :) + pars.Jta  * s.n.SPa(s.t, :);
    % Pitch-to-pitch 
    xPei = sum(pars.Jei * pars.Cei .* s.p.Sn(indPS) + ...
               pars.Jai * pars.Cei .* s.p.Sa(indPS), 2)';
    xPii = sum(pars.Jii * pars.Cii .* s.p.Sg(indPS), 2)';
    xPee = sum(pars.Jee * pars.Cee .* s.p.Sn(indPS) + ...
               pars.Jae * pars.Cee .* s.p.Sa(indPS), 2)';
    xPie = sum(pars.Jie * pars.Cie .* s.p.Sg(indPS), 2)';
    % Q-to-pitch inputs
    zPee = pars.Jqpen * s.q.Sn(indQS, :) + pars.Jqpea * s.q.Sa(indQS, :);
    zPei = pars.Jqpin * s.q.Sn(indQS, :) + pars.Jqpia * s.q.Sa(indQS, :);
    zPie = pars.Jqpeg * s.q.Sg(indQS, :);
    zPii = pars.Jqpig * s.q.Sg(indQS, :);

    % Q inputs 
    % Pitch-to-Q
    zQee = pars.Jpqen * s.p.Sn(indQS, :) + pars.Jpqea * s.p.Sa(indQS, :);
    zQei = pars.Jpqin * s.p.Sn(indQS, :) + pars.Jpqia * s.p.Sa(indQS, :);
    zQie = pars.Jpqeg * s.p.Sg(indQS, :);
    zQii = pars.Jpqig * s.p.Sg(indQS, :);
    % Q-to-Q (same AMPA conductivities as in pitch-to-pitch)
    xQei = sum(pars.Jqei * pars.Cqei .* s.q.Sn(indPS) + ...
               pars.Jqai * pars.Cqei .* s.q.Sa(indPS), 2)';
    xQii = sum(pars.Jqii * pars.Cqii .* s.q.Sg(indPS), 2)';
    xQee = sum(pars.Jqee * pars.Cqee .* s.q.Sn(indPS) + ...
               pars.Jqae * pars.Cqee .* s.q.Sa(indPS), 2)';
    xQie = sum(pars.Jqie * pars.Cqie .* s.q.Sg(indPS), 2)';

    % Population I0s
    extE = repmat(pars.I0e, size(xPee));
    extI = repmat(pars.I0i, size(xPii));

    % Net inputs
    xPi = extI + xPei - xPii + zPei - zPii - s.p.Ai(s.t, :);
    xPe = extE + xPee - xPie + zPee - zPie - s.p.Ae(s.t, :) + thalP;

    xQi = extI + xQei - xQii + zQei - zQii - s.q.Ai(s.t, :) + pars.QI0i;
    xQe = extE + xQee - xQie + zQee - zQie - s.q.Ae(s.t, :) + pars.QI0e;

end



function [dSn, dSp, dSq] = dSdt(s, r, pars)

    n = pars.sigma * randn([1, 13]);

    leakPI = -s.p.Sg(s.t, :)  / pars.tauGABA;
    leakPE = -s.p.Sn(s.t, :)  / pars.tauNMDA;
    leakPA = -s.p.Sa(s.t, :)  / pars.tauAMPA;

    dSp.g  = leakPI + (1 / 1000) * s.p.Hi(s.t, :) + n(4);
    dSp.n  = leakPE + pars.gamma * (1-s.p.Sn(s.t, :)).*s.p.He(s.t, :) + n(5);
    dSp.a  = leakPA + (1 / 1000) * s.p.He(s.t, :) + n(6);
    
    leakQI = -s.q.Sg(s.t, :) / pars.tauGABA;
    leakQE = -s.q.Sn(s.t, :) / pars.tauNMDA;
    leakQA = -s.q.Sa(s.t, :) / pars.tauAMPA;

    dSq.g  = leakQI + (1 / 1000) * s.q.Hi(s.t, :) + n(7);
    dSq.n  = leakQE + pars.gamma * (1-s.q.Sn(s.t, :)).*s.q.He(s.t, :) + n(8);
    dSq.a  = leakQA + (1 / 1000) * s.q.He(s.t, :) + n(9);

    leakNPn = -s.n.SPn(s.t, :) / pars.tauNMDA;
    leakNPa = -s.n.SPa(s.t, :) / pars.tauAMPA;

    dSn.Pn = leakNPn + pars.gamma * (1-s.n.SPn(s.t, :)).*r.A(s.t, :) + n(10);
    dSn.Pa = leakNPa + (1 / 1000) * r.A(s.t, :) + n(11);

end



function [dAPi, dAPe, dAQi, dAQe] = dAdt(s, pars)

	dAPi = -s.p.Ai(s.t, :) / pars.tauAdapt + pars.adaptStr * s.p.Hi(s.t, :);
    dAPe = -s.p.Ae(s.t, :) / pars.tauAdapt + pars.adaptStr * s.p.He(s.t, :);
    dAQi = -s.q.Ai(s.t, :) / pars.tauAdapt + pars.adaptStr * s.q.Hi(s.t, :);
    dAQe = -s.q.Ae(s.t, :) / pars.tauAdapt + pars.adaptStr * s.q.He(s.t, :);

end



function [dHPi, dHPe, dHQi, dHQe] = dHdt(s, pars)

    [phiIpe, tauPe] = phi(s.p.xe(s.t, :), s.p.He(s.t, :), 'e', pars);
    [phiIpi, tauPi] = phi(s.p.xi(s.t, :), s.p.Hi(s.t, :), 'i', pars);
    [phiIqe, tauQe] = phi(s.q.xe(s.t, :), s.q.He(s.t, :), 'e', pars);
    [phiIqi, tauQi] = phi(s.q.xi(s.t, :), s.q.Hi(s.t, :), 'i', pars);

    dHPe = (phiIpe - s.p.He(s.t, :)) ./ tauPe; 
    dHPi = (phiIpi - s.p.Hi(s.t, :)) ./ tauPi;
    dHQe = (phiIqe - s.q.He(s.t, :)) ./ tauQe; 
    dHQi = (phiIqi - s.q.Hi(s.t, :)) ./ tauQi;

end



function [phi, tau] = phi(x, H, typ, pars)

    Delta = 1; 
    a = pars.(['a', typ]);
    b = pars.(['b', typ]);
    d = pars.(['d', typ]);

    y = a * x - b;
    phi = y ./ (1 - exp(-d * y));
    
    if pars.tauEff
        dp = a * (1 ./ (y + eps) + d ./ (1 - exp(d * y))) .* phi; 
        tau = pars.tauPop * min(1, Delta * dp' ./ H')';
    else
        tau = pars.tauPop;
    end

    tau = max(1000 / pars.cortFs, tau); % tau cannot be smaller than dt

end
