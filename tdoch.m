function [s, r, lagSpace, timeSpace] = tdoch(pars, parsing)

    %if nargin == 0; pars = loadParameters(); end
    if nargin == 0; pars = loadParameters(); end
    if nargin < 2;  parsing = 0; end;

    % Sound --> auditory nerve spiking probabilities
    if pars.verb, fprintf('Computing thalamic input...\n'); end;

    r.lagSpace = binSpace(pars);
    r.freqSpace = 1000 ./ r.lagSpace; 

    if ischar(parsing)
        [timeSpace, r.A, r.n, r.b] = parseThalamic(parsing);
    else
        [timeSpace, r.A, r.n, r.b] = pyThalamic(r.lagSpace, pars);
    end

    % Subcortical libs are in seconds but cortical libs are in ms (sorry)
    r.timeSpace = timeSpace * 1000; 
    
    % Computing system evolution
    if not(pars.onlySubcort)
        s = tdochCortex(r, pars);
    else
        s = 0;
    end        

    timeSpace = r.timeSpace;
    lagSpace  = r.lagSpace;

end



function [timeSpace, A, n, b] = pyThalamic(lagSpace, pars)
 
    % parse filename is randomized to allow parallel computations
    chart = char(['A':'Z' 'a':'z' '0':'9']);
    parseID = ['pyparse' chart(ceil(length(chart) * rand(1, 4)))];

    save([parseID, 'In.mat'], 'pars', 'lagSpace');
 
    % MATLAB accessed bash didn't include many of my bash paths for some 
    % reason, so I had to use ipython to allow python access to the packages
    % installed through pip and add anaconda's libraries manually to allow 
    % bash to access to ipython. Adjust these lines as necessary!
    %python = '/home/tabs/Apps/anaconda/bin/ipython --colors=NoColor';
    python = 'python';

    system([python, ' subthalamic.py ', parseID]);

    [timeSpace, A, n, b] = parseThalamic([parseID, 'Out.mat']);
    
    delete([parseID, 'Out.mat']);

end



function [timeSpace, A, n, b] = parseThalamic(parsing)
 
    pyparse = load(parsing);
    timeSpace = pyparse.timeSpace;
    A = pyparse.A;    
    n = pyparse.n;
    b = pyparse.b;

end



function lagSpace = binSpace(pars)

    lagMin = 1000 / (pars.freqInterval(2));
    lagMax = 1000 / (pars.freqInterval(1));
    lagSpace = linspace(lagMax, lagMin, pars.N)';

end


