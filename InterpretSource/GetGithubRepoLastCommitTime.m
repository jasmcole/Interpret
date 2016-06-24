function date = GetGithubRepoLastCommitTime(repo)

optionsText = weboptions('username','jasmcole', 'password', 'dd20afd5570e01e418a0e59ea906b0e8e71f2652');
resp = webread(['https://api.github.com/repos/jasmcole/' repo '/git/refs/heads/master'], optionsText);
commit = webread(resp.object.url);
datetime = commit.author.date;
date = datetime(1:10);

end