function date = GetGithubRepoLastCommitTime(repo, srcPrefix)

token = load([srcPrefix 'update.mat']);
optionsText = weboptions('username','jasmcole', 'password', token.token);
resp = webread(['https://api.github.com/repos/jasmcole/' repo '/git/refs/heads/master'], optionsText);
commit = webread(resp.object.url);
datetime = commit.author.date;
date = datetime(1:10);

end