# 自动化连接 master 节点
ssh_from=$(who -m|awk -F "(" '{print $2}'|sed 's/)//g')
# 如果 from 为 server01 往其他节点跳转 则不执行
#if echo "$ssh_from" | grep -q "server"; then
#    return
#else  # 如果不是来源于 server 节点的连接 则跳转到 server
#    if [ "$(hostname)" != "server01" ]; then
#        user=***
#        password=***
#        /gpfs/oe-scrna/zhengfuxing/software/bin/sshpass -p $password ssh ${user}@10.100.100.1
#    fi
#fi
