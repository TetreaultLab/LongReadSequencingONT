# 1. Install the inotify tool for watch
sudo apt-get install inotify-tools
# 2. Add the watch_new_seq.sh to your home (~ or /home/username) and activate it with "sudo chmod +x ~/watch_new_seq.sh"
# 3. Add the config_prompt.sh to your desktop (~/Desktop/.) and activate it with "sudo chmod +x ~/Desktop/config_prompt.sh"
# 4. Add the service config to your user configs (meaning if you have multiple user sessions, needs to be added multiple times 
# For this just move seqwatch.service to ~/.config/systemd/user/ and activate it with "sudo chmod +x ~/.config/systemd/user/seqwatch.service"
# 5. Reboot (As the service triggers on restart) 
# 
# Important note : If modifications are made to the scripts, reboot may be required for the service to pick up the new scripts (restarting the deamon may not suffice) but can be tried:
systemctl --user daemon-reload
systemctl --user restart seqwatch.service
