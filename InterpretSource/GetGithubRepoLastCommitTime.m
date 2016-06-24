function date = GetGithubRepoLastCommitTime(repo)

optionsText = weboptions('username','jasmcole', 'password', 'f56003177f78f893519e7eabc0a77957d3c86bc4');
resp = webread(['https://api.github.com/repos/jasmcole/' repo '/git/refs/heads/master'], optionsText);
commit = webread(resp.object.url);
datetime = commit.author.date;
date = datetime(1:10);

end