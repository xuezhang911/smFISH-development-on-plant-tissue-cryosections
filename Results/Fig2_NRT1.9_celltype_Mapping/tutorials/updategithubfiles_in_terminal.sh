#initiate the local directory as a Git repository
git init -b main
#clone the github repository to local computal using git
git clone https://github.com/xuezhang911/tutorials.git
# enter the lcoal folder
cd tutorials
#update the github website
git pull origin main 
# put a new bash script in this folder
nano test.sh
ls
git status
# add the modified file to my repository
git add test.sh
git commit -a -m "initial commit"
git push origin main
#optionally git git add . && git commit -m "initial commit" git push orgin mian 
# update the github website and you will see the test.sh in tutorials
