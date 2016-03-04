function date = GetGithubRepoLastCommitTime(repo)

resp = webread(['https://api.github.com/repos/jasmcole/' repo '/git/refs/heads/master']);
commit = webread(resp.object.url);
datetime = commit.author.date;
date = datetime(1:10);

end