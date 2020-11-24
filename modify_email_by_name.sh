
#!/bin/sh
#
 
git filter-branch --force --env-filter '
    if [ "$GIT_COMMITTER_NAME" = "Yunlu Gao" ];
    then
        GIT_COMMITTER_NAME="B179895-2020";
        GIT_COMMITTER_EMAIL="<email_not_provided>";
        GIT_AUTHOR_NAME="B179895-2020";
        GIT_AUTHOR_EMAIL="<email_not_provided>";
    fi' -- --all
